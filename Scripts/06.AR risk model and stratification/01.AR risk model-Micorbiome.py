"""
Usage:
  python 01.AR risk model-Micorbiome.py
  python 01.AR risk model-Micorbiome.py -i Target_microbiome.csv -o .
  python 01.AR risk model-Micorbiome.py --input Target_microbiome.csv --outdir results --label-column Group
"""

import argparse
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import shap
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, classification_report, roc_auc_score, roc_curve, auc
from sklearn.model_selection import GridSearchCV, StratifiedKFold, train_test_split
from sklearn.preprocessing import LabelEncoder

matplotlib.rcParams["font.sans-serif"] = ["Arial"]


def parse_args():
    parser = argparse.ArgumentParser(description="Train a random forest microbiome AR classifier and export ROC/SHAP outputs.")
    parser.add_argument("-i", "--input", default="Target_microbiome.csv", help="Input CSV file (default: Target_microbiome.csv)")
    parser.add_argument("-o", "--outdir", default=".", help="Output directory (default: current directory)")
    parser.add_argument("-l", "--label-column", default="Group", help="Label column name (default: Group)")
    return parser.parse_args()


args = parse_args()
input_path = Path(args.input)
outdir = Path(args.outdir)

if not input_path.is_file():
    raise FileNotFoundError(f"Input file not found: {input_path}")

outdir.mkdir(parents=True, exist_ok=True)

# 1. Import the dataset
df = pd.read_csv(input_path)

if args.label_column not in df.columns:
    raise ValueError(f"Label column not found: {args.label_column}")

feature_cols = df.columns[2:]
if len(feature_cols) == 0:
    raise ValueError("No feature columns found after the first two columns.")

# 2. Split X and y
X = df.loc[:, feature_cols]
y = df[args.label_column]

# 3. Encode group labels (HC=0, AR=1)
y_encoded = y.astype(pd.CategoricalDtype(categories=["HC", "AR"], ordered=True)).cat.codes
if (y_encoded < 0).any():
    raise ValueError("Label encoding produced NA values. Ensure the label column only contains HC and AR.")
if y_encoded.nunique() < 2:
    raise ValueError("At least two groups are required in the label column.")

# 4. Model training

# Split dataset
X_train, X_test, y_train, y_test = train_test_split(X, y_encoded, test_size=0.2, random_state=123, stratify=y_encoded)

# Setting the hyperparameter search space
param_grid = {
    "n_estimators": [100, 200, 300],
    "max_depth": [5, 10, 15, None],
    "min_samples_leaf": [1, 2, 4],
    "max_features": ["sqrt", "log2"]
}

# Optimize hyperparameters using GridSearchCV (5-fold cross-validation)
rf = RandomForestClassifier(random_state=123)
grid_search = GridSearchCV(
    rf,
    param_grid,
    cv=5,
    scoring="roc_auc",
    n_jobs=-1,
    verbose=1
)
grid_search.fit(X_train, y_train)

# Optimal Model
best_model = grid_search.best_estimator_
print("✅ Best Parameters:", grid_search.best_params_)

# Model Evaluation
y_pred = best_model.predict(X_test)
y_prob = best_model.predict_proba(X_test)[:, 1]

acc = accuracy_score(y_test, y_pred)
auc = roc_auc_score(y_test, y_prob)
report = classification_report(y_test, y_pred, target_names=LabelEncoder().fit(y).classes_)

print(f"\n✅ Test Accuracy: {acc:.4f}")
print(f"✅ Test AUC: {auc:.4f}")
print("✅ Classification Report:\n", report)

# 5. ROC Curve with 5-fold Cross-Validation
cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=123)

tprs = []
aucs = []
mean_fpr = np.linspace(0, 1, 100)

fig, ax = plt.subplots(figsize=(6, 6))

for i, (train, test) in enumerate(cv.split(X, y_encoded)):
    model = best_model.fit(X.iloc[train], y_encoded.iloc[train])
    y_prob = model.predict_proba(X.iloc[test])[:, 1]

    fpr, tpr, _ = roc_curve(y_encoded.iloc[test], y_prob)
    roc_auc = auc(fpr, tpr)
    aucs.append(roc_auc)
    ax.plot(fpr, tpr, lw=1, alpha=0.3, label=f"Fold {i+1} (AUC = {roc_auc:.2f})")

    tprs.append(np.interp(mean_fpr, fpr, tpr))
    tprs[-1][0] = 0.0

mean_tpr = np.mean(tprs, axis=0)
mean_tpr[-1] = 1.0
mean_auc = auc(mean_fpr, mean_tpr)
std_auc = np.std(aucs)

ax.plot(
    mean_fpr,
    mean_tpr,
    color="b",
    label=r"Mean ROC (AUC = %0.2f $\pm$ %0.2f)" % (mean_auc, std_auc),
    lw=2,
    alpha=0.8,
)

std_tpr = np.std(tprs, axis=0)
tpr_upper = np.minimum(mean_tpr + std_tpr, 1)
tpr_lower = np.maximum(mean_tpr - std_tpr, 0)
ax.fill_between(
    mean_fpr,
    tpr_lower,
    tpr_upper,
    color="grey",
    alpha=0.2,
    label=r"$\pm$ 1 std. dev.",
)

ax.plot([0, 1], [0, 1], linestyle="--", color="r", lw=2)
ax.set_xlim([0, 1])
ax.set_ylim([0, 1.05])
ax.set_xlabel("False Positive Rate")
ax.set_ylabel("True Positive Rate")
ax.set_title("5-Fold Cross-Validated ROC (Best RF Model)")
ax.legend(loc="lower right")
plt.tight_layout()
plt.savefig(outdir / "ROC_curve_microbiome.pdf", format="pdf", bbox_inches="tight")
plt.show()
plt.close(fig)

# 6. SHAP Analysis
explainer = shap.TreeExplainer(best_model, model_output="raw")
shap_values = explainer.shap_values(X)

mean_shap_vals = np.abs(shap_values[1]).mean(axis=0)
feature_importance = pd.DataFrame(
    {
        "feature": X_test.columns,
        "mean_abs_shap": mean_shap_vals,
    }
)
feature_importance.sort_values("mean_abs_shap", ascending=False, inplace=True)

feature_importance.to_csv(outdir / "Shap_feature_importance_microbiome.csv", index=False)
