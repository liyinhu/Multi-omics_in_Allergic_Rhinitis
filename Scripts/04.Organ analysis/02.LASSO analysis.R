# Usage:
#   Rscript 02.LASSO analysis.R
#   Rscript 02.LASSO analysis.R --expr expr.csv --organ-map organ_map.csv --organ-index organ_index.csv
#
# If no arguments are provided, the script falls back to interactive file selection.

library(glmnet)
library(caret)
library(dplyr)
library(ggplot2)
library(stringr)
library(tidyr)

reorder_within <- function(x, by, within, fun = mean, sep = "___", ...) {
  new_x <- paste(x, within, sep = sep)
  stats::reorder(new_x, by, FUN = fun, ...)
}

scale_y_reordered <- function(..., sep = "___") {
  scale_y_discrete(labels = function(x) sub(paste0(sep, ".+$"), "", x), ...)
}

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

cli_args <- commandArgs(trailingOnly = TRUE)
expr_path <- parse_cli_value(cli_args, c("--expr", "-e"))
organ_map_path <- parse_cli_value(cli_args, c("--organ-map", "-m"))
organ_index_path <- parse_cli_value(cli_args, c("--organ-index", "-i"))

if (is.null(expr_path)) {
  expr_path <- file.choose()
}
if (is.null(organ_map_path)) {
  organ_map_path <- file.choose()
}
if (is.null(organ_index_path)) {
  organ_index_path <- file.choose()
}

if (!file.exists(expr_path)) {
  stop(sprintf("Expression matrix not found: %s", expr_path))
}
if (!file.exists(organ_map_path)) {
  stop(sprintf("Protein-organ mapping file not found: %s", organ_map_path))
}
if (!file.exists(organ_index_path)) {
  stop(sprintf("Organ index file not found: %s", organ_index_path))
}

# Import proteome matrix
expr_df <- read.csv(expr_path, header = TRUE, row.names = 1)

# Import protein-organ mapping matrix
organ_df <- read.csv(organ_map_path, header = TRUE)

# Import organ matrix
organ_index_df <- read.csv(organ_index_path, header = TRUE, row.names = 1)

required_organ_cols <- c("Tissue2", "Gene_ID")
missing_organ_cols <- setdiff(required_organ_cols, colnames(organ_df))
if (length(missing_organ_cols) > 0) {
  stop(sprintf(
    "Protein-organ mapping file is missing required columns: %s",
    paste(missing_organ_cols, collapse = ", ")
  ))
}

top_n <- 10
all_feature_results <- list()

# Processing
for (organ in colnames(organ_index_df)) {
  message("Organ: ", organ)

  # Isolate proteins
  protein_features <- organ_df %>%
    filter(.data$Tissue2 == organ) %>%
    pull(Gene_ID) %>%
    unique()
  protein_features <- intersect(protein_features, colnames(expr_df))

  if (length(protein_features) < 5) next

  # Construct the feature matrix
  X <- expr_df[, protein_features, drop = FALSE]
  y <- organ_index_df[[organ]]

  # Training/Test
  set.seed(123)
  train_idx <- createDataPartition(y, p = 0.7, list = FALSE)
  X_train <- X[train_idx, ]
  y_train <- y[train_idx]
  X_test <- X[-train_idx, ]
  y_test <- y[-train_idx]

  # Lasso model + fine tuning
  cvfit <- cv.glmnet(as.matrix(X_train), as.numeric(y_train), alpha = 1, nfolds = 5)
  best_lambda <- cvfit$lambda.min
  final_model <- glmnet(as.matrix(X), as.numeric(y), alpha = 1, lambda = best_lambda)

  # Get LASSO coeffecients
  coef_df <- coef(final_model)
  coef_table <- data.frame(
    feature = rownames(coef_df),
    coefficient = as.vector(coef_df)
  ) %>%
    filter(feature != "(Intercept)", coefficient != 0)

  top_features <- coef_table %>%
    mutate(abs_coef = abs(coefficient)) %>%
    arrange(desc(abs_coef)) %>%
    slice_head(n = top_n) %>%
    mutate(organ = organ) %>%
    select(organ, feature, coefficient)

  all_feature_results[[organ]] <- top_features
}

if (length(all_feature_results) == 0) {
  stop("No organ passed the feature threshold (length(protein_features) < 5).")
}

# Export results
final_df <- do.call(rbind, all_feature_results)
write.csv(final_df, "Top_Lasso_Features_by_Organ.csv", row.names = FALSE)

# Plotting
df <- final_df %>% mutate(feature = reorder_within(feature, coefficient, organ))

ggplot(df, aes(x = coefficient, y = feature, fill = organ, color = organ)) +
  geom_segment(aes(y = feature, yend = feature, x = 0, xend = coefficient), linewidth = 0.8, color = "grey80") +
  geom_point(size = 2) +
  facet_wrap(~organ, scales = "free") +
  scale_y_reordered() +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_fill_manual(values = c(
    "Organism" = "#D86C76",
    "ReproductiveF" = "#0D6CA6",
    "Hepatobiliary" = "#E6DB49",
    "Gastrointestinal" = "#CB5D17",
    "Brain" = "#14976F",
    "Endocrine" = "#58ABDB",
    "Lymphoidtissue" = "#C2749E",
    "ReproductiveM" = "#C36463",
    "Bonemarrow" = "#ABCCDE",
    "Urinary" = "#69BE92",
    "Pancreas" = "#DD9B15",
    "Adiposetissue" = "#F5CADE",
    "Skin" = "#ADD16A",
    "Proximal digestive tract" = "#F4AA8E",
    "Lung" = "#FBE3D8",
    "Muscle" = "#BBB7D6",
    "Eye" = "#8ACBC0"
  )) +
  scale_color_manual(values = c(
    "Organism" = "#D86C76",
    "ReproductiveF" = "#0D6CA6",
    "Hepatobiliary" = "#E6DB49",
    "Gastrointestinal" = "#CB5D17",
    "Brain" = "#14976F",
    "Endocrine" = "#58ABDB",
    "Lymphoidtissue" = "#C2749E",
    "ReproductiveM" = "#C36463",
    "Bonemarrow" = "#ABCCDE",
    "Urinary" = "#69BE92",
    "Pancreas" = "#DD9B15",
    "Adiposetissue" = "#F5CADE",
    "Skin" = "#ADD16A",
    "Proximal digestive tract" = "#F4AA8E",
    "Lung" = "#FBE3D8",
    "Muscle" = "#BBB7D6",
    "Eye" = "#8ACBC0"
  ))
