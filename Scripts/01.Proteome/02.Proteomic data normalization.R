# Usage:
#   Rscript "02.Proteomic data normalization.R"
#   Rscript "02.Proteomic data normalization.R" --protein protein_matrix.txt --info sample_info.txt --outdir results

library(plotly)
library(robustbase)
library(vegan)

parse_cli_value <- function(args, flags) {
  idx <- match(flags, args)
  idx <- idx[!is.na(idx)][1]
  if (is.na(idx)) {
    return(NULL)
  }
  if (idx == length(args)) {
    stop(sprintf("Missing value after %s", args[[idx]]))
  }
  args[[idx + 1]]
}

read_delim_checked <- function(path, label) {
  if (!file.exists(path)) {
    stop(sprintf("%s file not found: %s", label, path))
  }
  read.delim(path, header = TRUE, row.names = 1, check.names = FALSE)
}

args <- commandArgs(trailingOnly = TRUE)
protein_path <- parse_cli_value(args, c("--protein", "-p"))
info_path <- parse_cli_value(args, c("--info", "-i"))
outdir <- parse_cli_value(args, c("--outdir", "-o"))

if (is.null(outdir) || !nzchar(outdir)) {
  outdir <- "."
}
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Input proteomic and background files
if (is.null(protein_path) || !nzchar(protein_path)) {
  protein_path <- file.choose()
}
if (is.null(info_path) || !nzchar(info_path)) {
  info_path <- file.choose()
}

data <- read_delim_checked(protein_path, "Protein matrix")
info <- read_delim_checked(info_path, "Sample info")

if (nrow(data) == 0 || ncol(data) == 0) {
  stop("Protein matrix must contain samples in rows and proteins in columns.")
}
if (nrow(info) == 0 || ncol(info) == 0) {
  stop("Sample info table is empty.")
}
if (is.null(rownames(data)) || is.null(rownames(info))) {
  stop("Protein matrix and sample info must both have sample IDs in the first column / row names.")
}

required_info_cols <- c("Group", "Gender", "Age", "Year")
missing_info_cols <- setdiff(required_info_cols, colnames(info))
if (length(missing_info_cols) > 0) {
  stop(sprintf(
    "Sample info table is missing required columns: %s",
    paste(missing_info_cols, collapse = ", ")
  ))
}

common_samples <- intersect(rownames(data), rownames(info))
if (length(common_samples) == 0) {
  stop("Protein matrix and sample info have no overlapping sample IDs.")
}

data <- data[common_samples, , drop = FALSE]
info <- info[common_samples, , drop = FALSE]

if (!identical(rownames(data), rownames(info))) {
  stop("Protein matrix and sample info could not be aligned to the same sample order.")
}

protein_numeric <- suppressWarnings(
  as.data.frame(lapply(data, function(x) as.numeric(as.character(x))), check.names = FALSE)
)
rownames(protein_numeric) <- rownames(data)

bad_protein_cols <- names(data)[vapply(
  seq_along(data),
  function(i) any(is.na(protein_numeric[[i]]) & !is.na(data[[i]])),
  logical(1)
)]
if (length(bad_protein_cols) > 0) {
  stop(sprintf(
    "Protein matrix must contain only numeric protein expression values. Invalid columns: %s",
    paste(bad_protein_cols, collapse = ", ")
  ))
}

data <- protein_numeric

age_raw <- info$Age
info$Age <- suppressWarnings(as.numeric(as.character(age_raw)))
if (any(is.na(info$Age) & !is.na(age_raw))) {
  stop("Age column must be numeric or convertible to numeric.")
}
if (all(is.na(info$Age))) {
  stop("Age column could not be converted to numeric.")
}

gender_raw <- info$Gender
gender_numeric <- suppressWarnings(as.numeric(as.character(gender_raw)))
gender_is_numeric <- !any(is.na(gender_numeric) & !is.na(gender_raw))

if (gender_is_numeric) {
  info$Gender <- gender_numeric
} else {
  info$Gender <- as.factor(gender_raw)
}

protein_count <- ncol(data)
if (protein_count != 362) {
  warning(sprintf(
    "Expected 362 protein columns based on the original script, but detected %d. Continuing with the detected protein columns.",
    protein_count
  ))
}

# Permutational multivariate analysis of variance (Before normalization)
permonova <- adonis(data ~ Group + Gender + Age + Year, data = info, permutations = 9999, method = "euclidean")

# Construct a table for correction coefficients and adjusted results
Adjusted_data <- data.frame(matrix(ncol = protein_count, nrow = nrow(data)))
colnames(Adjusted_data) <- colnames(data)
rownames(Adjusted_data) <- rownames(data)

Age_Sex_Effects_on_proteins <- data.frame(matrix(ncol = 8, nrow = protein_count))
colnames(Age_Sex_Effects_on_proteins) <- c(
  "Age_Estimate", "Age_SE", "Age_T_value", "Age_P_value",
  "Sex_Estimate", "Sex_SE", "Sex_T_value", "Sex_P_value"
)
rownames(Age_Sex_Effects_on_proteins) <- colnames(data)

# Data normalization: linear regression to correct the effects of age and sex
Merge_data <- cbind(data, info)

for (k in seq_len(protein_count)) {
  protein_name <- colnames(data)[k]
  formula_text <- sprintf("`%s` ~ Age + Gender", protein_name)
  Protein_age_sex_test <- lmrob(stats::as.formula(formula_text), data = Merge_data, k.max = 900000)

  coef_summary <- summary(Protein_age_sex_test)$coefficients
  age_row <- which(rownames(coef_summary) == "Age")
  gender_rows <- grep("^Gender", rownames(coef_summary))

  if (length(age_row) == 1) {
    Age_Sex_Effects_on_proteins[k, 1:4] <- as.vector(coef_summary[age_row, 1:4])
  }
  if (length(gender_rows) >= 1) {
    Age_Sex_Effects_on_proteins[k, 5:8] <- as.vector(coef_summary[gender_rows[1], 1:4])
    if (length(gender_rows) > 1) {
      warning(sprintf(
        "Protein %s has multiple Gender coefficients; only the first Gender coefficient is written to AgeSex_Effects.csv.",
        protein_name
      ))
    }
  }

  model_matrix <- model.matrix(Protein_age_sex_test)
  coefficient_vector <- stats::coef(Protein_age_sex_test)
  non_intercept_idx <- which(colnames(model_matrix) != "(Intercept)")
  adjustment <- model_matrix[, non_intercept_idx, drop = FALSE] %*% coefficient_vector[non_intercept_idx]
  Adjusted_data[, k] <- Merge_data[[protein_name]] - as.numeric(adjustment)
}

# Export correction coefficients and adjusted results
write.csv(Adjusted_data, file.path(outdir, "Agesex_adjusted_protein_levels.csv"), row.names = TRUE)
write.csv(Age_Sex_Effects_on_proteins, file.path(outdir, "AgeSex_Effects.csv"), row.names = TRUE)

# Permutational multivariate analysis of variance (After normalization)
permonova2 <- adonis(Adjusted_data ~ Group + Gender + Age + Year, data = info, permutations = 9999, method = "euclidean")
