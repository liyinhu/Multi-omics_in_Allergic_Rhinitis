# Usage:
#   python organ_risk_index.py -p Proteome/Adjusted/All.proteome.res_adjusted.csv -o Protein_Organ.v2.list -i Info.csv -r Organ
#
# Required arguments:
#   -p / --protein     Protein expression matrix (CSV, rows = samples, columns = proteins)
#   -o / --organ       Protein-to-organ mapping table (TSV, must contain Organ and Protein columns)
#   -i / --info        Sample metadata table (CSV, must contain Group column)
#   -r / --result      Output directory, created if it does not exist
#
# Optional argument:
#   -m / --metabolite   Kept for interface compatibility; currently not used by the script.

import argparse
import os
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from sklearn.calibration import CalibratedClassifierCV
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import GridSearchCV
from sklearn.utils import resample
from tqdm import tqdm


def build_parser():
    parser = argparse.ArgumentParser(
        description="Estimate organ-level risk indices from protein expression data.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  python organ_risk_index.py -p Proteome/Adjusted/All.proteome.res_adjusted.csv "
            "-o Protein_Organ.v2.list -i Info.csv -r Organ\n"
            "  python organ_risk_index.py --protein protein.csv --organ Protein_Organ.v2.list "
            "--info Info.csv --result Organ"
        ),
    )
    parser.add_argument("-p", "--protein", dest="protein", action="store", type=str, required=True, help="protein file")
    parser.add_argument("-m", "--metabolite", dest="metabolite", action="store", type=str, required=False, help="metabolite file")
    parser.add_argument("-o", "--organ", dest="organ", action="store", type=str, required=True, help="protein to organ file")
    parser.add_argument("-i", "--info", dest="info", action="store", type=str, required=True, help="sample info")
    parser.add_argument("-r", "--result", dest="result", action="store", type=str, required=False, default="./", help="result directory")
    return parser


def validate_file(path_str, label):
    path = Path(path_str)
    if not path.is_file():
        raise FileNotFoundError(f"{label} not found: {path}")
    return path


def main():
    args = build_parser().parse_args()

    protein_path = validate_file(args.protein, "Protein file")
    organ_path = validate_file(args.organ, "Protein-to-organ file")
    info_path = validate_file(args.info, "Sample info file")

    if args.metabolite:
        validate_file(args.metabolite, "Metabolite file")

    result_dir = Path(args.result)
    result_dir.mkdir(parents=True, exist_ok=True)

    # protein_expr_df, protein2organ_df, sample_meta_df
    protein_expr_df = pd.read_csv(protein_path, sep=",", index_col=0, header=0, low_memory=False)  # Row: sample; Column: protein
    protein2organ_df = pd.read_csv(organ_path, sep="\t")
    sample_meta_df = pd.read_csv(info_path, sep=",", index_col=0, header=0, low_memory=False)  # Row: sample; Column: Age, Sex, Group, etc

    required_organ_cols = {"Organ", "Protein"}
    missing_organ_cols = required_organ_cols.difference(protein2organ_df.columns)
    if missing_organ_cols:
        raise ValueError(f"Protein-to-organ file is missing required columns: {sorted(missing_organ_cols)}")

    if "Group" not in sample_meta_df.columns:
        raise ValueError("Sample info file must contain a 'Group' column.")

    all_lower = all(isinstance(col, str) and col.islower() for col in protein_expr_df.columns)
    if all_lower:
        protein_expr_df.columns = protein_expr_df.columns.str.upper()
        print("Detected lowercase protein names; converted columns to uppercase.")

    # 基础准备
    n_bootstrap = 100
    n_jobs = 50
    results_mean_dict = {}
    results_se_dict = {}
    auc_mean_dict = {}
    auc_se_dict = {}

    # Isolate all organs
    organs = protein2organ_df["Organ"].dropna().unique()

    # ===================================
    # Parallel processing
    # ===================================
    def process_organ(organ):
        organ_proteins = protein2organ_df.query("Organ == @organ")["Protein"]
        organ_proteins = [p for p in organ_proteins if p in protein_expr_df.columns]
        if len(organ_proteins) < 2:
            return organ, None, None, None, None

        X_full = protein_expr_df[organ_proteins]
        y_full = sample_meta_df.loc[X_full.index, "Group"].map({"2AR": 1, "1HC": 0})
        if y_full.isna().any():
            return organ, None, None, None, None
        if y_full.value_counts().min() < 2:
            return organ, None, None, None, None

        pred_list = []
        auc_list = []
        for i in range(n_bootstrap):
            try:
                X_boot, y_boot = resample(X_full, y_full, replace=True, n_samples=len(X_full), random_state=i)
                param_grid = {"n_estimators": [100], "max_depth": [5], "min_samples_split": [2]}
                rf = RandomForestClassifier(random_state=i)
                rf_grid = GridSearchCV(rf, param_grid, cv=2, scoring="roc_auc", n_jobs=1)
                rf_grid.fit(X_boot, y_boot)
                calibrated_rf = CalibratedClassifierCV(rf_grid.best_estimator_, method="sigmoid", cv=2)
                calibrated_rf.fit(X_boot, y_boot)
                probs = calibrated_rf.predict_proba(X_full)[:, 1]
                pred_list.append(probs)

                auc = roc_auc_score(y_full, probs)
                auc_list.append(auc)
            except Exception:
                continue

        if len(pred_list) == 0:
            return organ, None, None, None, None

        pred_array = np.vstack(pred_list)
        return organ, pred_array.mean(axis=0), pred_array.std(axis=0), np.mean(auc_list), np.std(auc_list)

    # ===================================
    # Process
    # ===================================
    results = Parallel(n_jobs=n_jobs)(delayed(process_organ)(organ) for organ in tqdm(organs))

    # Get results
    for organ, mean_preds, std_preds, auc_mean, auc_se in results:
        if mean_preds is not None:
            results_mean_dict[organ] = mean_preds
            results_se_dict[organ] = std_preds
            auc_mean_dict[organ] = auc_mean
            auc_se_dict[organ] = auc_se

    results_mean = pd.DataFrame(results_mean_dict, index=protein_expr_df.index)
    results_se = pd.DataFrame(results_se_dict, index=protein_expr_df.index)
    results_mean.to_csv(result_dir / f"organ_risk_index_mean.{n_bootstrap}.csv")
    results_se.to_csv(result_dir / f"organ_risk_index_se.{n_bootstrap}.csv")

    auc_df = pd.DataFrame({
        "Organ": list(auc_mean_dict.keys()),
        "AUC_Mean": list(auc_mean_dict.values()),
        "AUC_SE": list(auc_se_dict.values())
    })
    auc_df.to_csv(result_dir / f"organ_risk_index_auc.{n_bootstrap}.csv", index=False)


if __name__ == "__main__":
    main()
