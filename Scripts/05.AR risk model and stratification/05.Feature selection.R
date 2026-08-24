# Usage:
#   Rscript "05.Feature selection.R"
#   Rscript "05.Feature selection.R" --input Target_multi_omics.csv --outdir results --label-column Group
#   Rscript "05.Feature selection.R" --input Target_multi_omics.csv --outdir results --label-column V14

library(glmnet)
library(Metrics)
library(ModelMetrics)
library(caret)

get_arg_value <- function(args, keys) {
  for (i in seq_along(args)) {
    if (args[[i]] %in% keys && i < length(args)) {
      return(args[[i + 1]])
    }
    for (key in keys) {
      if (startsWith(args[[i]], paste0(key, "="))) {
        return(sub(paste0("^", key, "="), "", args[[i]]))
      }
    }
  }
  return(NULL)
}

args <- commandArgs(trailingOnly = TRUE)
input_file <- get_arg_value(args, c("--input", "-i"))
outdir <- get_arg_value(args, c("--outdir", "-o"))
label_column <- get_arg_value(args, c("--label-column", "-l"))

if (is.null(outdir) || !nzchar(outdir)) {
  outdir <- "."
}
if (is.null(label_column) || !nzchar(label_column)) {
  label_column <- "Group"
}

if (is.null(input_file) || !nzchar(input_file)) {
  input_file <- file.choose()
}

if (!file.exists(input_file)) {
  stop(sprintf("Input file not found: %s", input_file))
}

if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
}

# Import multi-omics data (121 features with group information)
data <- read.delim(input_file, header = TRUE, row.names = 1, check.names = FALSE)

label_col_index <- NULL
if (label_column %in% names(data)) {
  label_col_index <- which(names(data) == label_column)[1]
} else if (ncol(data) >= 14) {
  label_col_index <- 14
  message(sprintf("Label column '%s' not found; falling back to column 14 ('%s').", label_column, names(data)[14]))
} else {
  stop(sprintf("Label column '%s' not found and data has fewer than 14 columns.", label_column))
}

if (ncol(data) <= 1) {
  stop("No feature columns found.")
}

feature_cols <- setdiff(seq_len(ncol(data)), label_col_index)
if (length(feature_cols) == 0) {
  stop("No feature columns found after excluding the label column.")
}

# Training and Test dataset
set.seed(123)
group_vec <- data[[label_col_index]]
if (length(unique(na.omit(group_vec))) < 2) {
  stop("At least two groups are required in the label column.")
}

partition_index <- createDataPartition(group_vec, p = 0.7, list = FALSE)
train_d <- data[partition_index, , drop = FALSE]
test_d <- data[-partition_index, , drop = FALSE]

# Training LASSO Model
lambdas <- seq(0, 2, length.out = 100)
X_raw <- train_d[, feature_cols, drop = FALSE]
feature_numeric <- suppressWarnings(data.frame(lapply(X_raw, function(x) as.numeric(as.character(x))), check.names = FALSE))

non_numeric_cols <- names(X_raw)[vapply(seq_along(X_raw), function(i) any(is.na(feature_numeric[[i]]) & !is.na(X_raw[[i]])), logical(1))]
if (length(non_numeric_cols) > 0) {
  stop(sprintf("The following feature columns cannot be converted to numeric: %s", paste(non_numeric_cols, collapse = ", ")))
}

if (ncol(feature_numeric) == 0) {
  stop("No feature columns found.")
}

X <- as.matrix(feature_numeric)
Y <- train_d[[label_col_index]]
if (length(unique(na.omit(Y))) < 2) {
  stop("At least two groups are required after splitting the training set.")
}

set.seed(123)
lasso_model <- cv.glmnet(X, Y, alpha = 1, lambda = lambdas, nfolds = 5)

# Plotting
pdf(file.path(outdir, "Lasso_feature_selection.pdf"), width = 8, height = 10)
plot(lasso_model)
plot(lasso_model$glmnet.fit, "lambda", label = TRUE)
dev.off()

# Select the optimal regularization parameters, train the optimal model, and output the Lasso coefficients of each variable.
lasso_min <- lasso_model$lambda.min
lasso_best <- glmnet(X, Y, alpha = 1, lambda = lasso_min)
write.csv(as.matrix(coef(lasso_best)), file.path(outdir, "Lasso_feature.csv"))
