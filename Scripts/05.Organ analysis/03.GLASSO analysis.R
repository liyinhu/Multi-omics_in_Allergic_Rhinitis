# Usage:
#   Rscript 03.GLASSO analysis.R
#   Rscript 03.GLASSO analysis.R --organ organ.csv --info info.tsv --outdir results
#
# If no arguments are provided, the script falls back to interactive file selection.

library(qgraph)
library(bootnet)
library(dplyr)
library(tidyr)
library(tidyverse)
library(corrplot)
library(RColorBrewer)

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
organ_path <- parse_cli_value(cli_args, c("--organ", "-o"))
info_path <- parse_cli_value(cli_args, c("--info", "-i"))
outdir <- parse_cli_value(cli_args, c("--outdir", "-r"))
if (is.null(outdir)) {
  outdir <- "."
}

if (is.null(organ_path)) {
  organ_path <- file.choose()
}
if (is.null(info_path)) {
  info_path <- file.choose()
}

if (!file.exists(organ_path)) {
  stop(sprintf("Organ matrix not found: %s", organ_path))
}
if (!file.exists(info_path)) {
  stop(sprintf("Sample info file not found: %s", info_path))
}

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

## 1. Import organ matrix
organ <- read.csv(organ_path, row.names = 1, check.names = FALSE)

# Import group information
info <- read.delim(info_path, header = TRUE, check.names = FALSE)  # Sample | Group (AR / HC)

required_info_cols <- c("ID", "Group")
missing_info_cols <- setdiff(required_info_cols, colnames(info))
if (length(missing_info_cols) > 0) {
  stop(sprintf(
    "Sample info file is missing required columns: %s",
    paste(missing_info_cols, collapse = ", ")
  ))
}

# Select common samples in organ-file order
common_samples <- rownames(organ)[rownames(organ) %in% info$ID]
if (length(common_samples) == 0) {
  stop("No overlapping sample IDs were found between the organ matrix and sample info.")
}
organ <- organ[common_samples, , drop = FALSE]
info  <- info[match(common_samples, info$ID), , drop = FALSE]

group <- info$Group

## 2. Z-score normalization
organ_z <- scale(organ)
organ_z <- as.data.frame(organ_z)

## 3. Grouping
organ_AR <- organ_z[group == "2AR", , drop = FALSE]
organ_HC <- organ_z[group == "1HC", , drop = FALSE]

if (nrow(organ_AR) < 2 || nrow(organ_HC) < 2) {
  stop("Each group must contain at least 2 samples after alignment for bootnet to run.")
}

if (all(apply(organ_AR, 2, sd, na.rm = TRUE) == 0) || all(apply(organ_HC, 2, sd, na.rm = TRUE) == 0)) {
  stop("One of the groups has zero-variance features only; bootnet cannot proceed.")
}

## 4. Bootstrap
## 4.1 Network estimation
estimate_glasso <- function(data) {
  EBICglasso(cor(data), n = nrow(data), gamma = 0.5)
}

## 4.2 Bootstrap
set.seed(123)

boot_AR <- bootnet(
  organ_AR,
  default = "EBICglasso",
  nBoots = 1000,
  type = "nonparametric",
  tuning = 0.75
)

boot_HC <- bootnet(
  organ_HC,
  default = "EBICglasso",
  nBoots = 1000,
  type = "nonparametric",
  tuning = 0.75
)

## 5. Isolate bootstrap results
## 5.1 Isolate bootTable for AR and HC
boot_AR_edges <- boot_AR$bootTable %>%
  filter(type == "edge") %>%
  rename(value_AR = value)

boot_HC_edges <- boot_HC$bootTable %>%
  filter(type == "edge") %>%
  rename(value_HC = value)

if (nrow(boot_AR_edges) == 0 || nrow(boot_HC_edges) == 0) {
  stop("No bootstrap edge results were returned by bootnet.")
}

## 5.2 Summarize the bootstrap statistics for each edge
edge_AR <- boot_AR_edges %>%
  group_by(node1, node2) %>%
  summarise(
    mean = mean(value_AR, na.rm = TRUE),
    sd   = sd(value_AR, na.rm = TRUE),
    CI_lower = quantile(value_AR, 0.025, na.rm = TRUE),
    CI_upper = quantile(value_AR, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

edge_HC <- boot_HC_edges %>%
  group_by(node1, node2) %>%
  summarise(
    mean = mean(value_HC, na.rm = TRUE),
    sd   = sd(value_HC, na.rm = TRUE),
    CI_lower = quantile(value_HC, 0.025, na.rm = TRUE),
    CI_upper = quantile(value_HC, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

## 6. Difference network between AR and HC
## 6.1 Merge the bootstrap edge tables of AR and HC
boot_AR_edges <- boot_AR$bootTable %>%
  filter(type == "edge") %>%
  rename(value_AR = value)

boot_HC_edges <- boot_HC$bootTable %>%
  filter(type == "edge") %>%
  rename(value_HC = value)

## 6.2 Paired by edge and bootstrap
boot_diff <- inner_join(
  boot_AR_edges,
  boot_HC_edges,
  by = c("node1", "node2", "name")
)

## 6.3 Calculate bootstrap differences
boot_diff <- boot_diff %>%
  mutate(diff = value_AR - value_HC)

## 6.4 Calculate the bootstrap p-value for each edge
edge_diff <- boot_diff %>%
  group_by(node1, node2) %>%
  summarise(
    mean_diff = mean(diff),
    p = {
      prop_neg <- mean(diff <= 0)
      prop_pos <- mean(diff >= 0)
      p_raw <- 2 * min(prop_neg, prop_pos)
      min(p_raw, 1)
    },
    CI_lower = quantile(diff, 0.025),
    CI_upper = quantile(diff, 0.975),
    .groups = "drop"
  )

if (nrow(edge_diff) == 0) {
  message("No differential edges were produced by the bootstrap comparison; writing empty outputs and skipping downstream plotting.")
}

## 6.5 Adjustement for multiple test
edge_diff$padj <- p.adjust(edge_diff$p, method = "BH")

## 7. Export results
write.csv(edge_AR, file.path(outdir, "glasso_edges_AR_bootstrap.csv"), row.names = FALSE)
write.csv(edge_HC, file.path(outdir, "glasso_edges_HC_bootstrap.csv"), row.names = FALSE)
write.csv(edge_diff, file.path(outdir, "glasso_edges_AR_vs_HC_diff.csv"), row.names = FALSE)

if (nrow(edge_diff) == 0) {
  empty_mat <- matrix(numeric(0), nrow = 0, ncol = 0)
  write.csv(empty_mat, file.path(outdir, "diff_edges_mat.csv"))
  write.csv(empty_mat, file.path(outdir, "diff_edges_p.csv"))
  write.csv(empty_mat, file.path(outdir, "diff_edges_padj.csv"))
} else {
## 8. Convert edge_diff into a matrix
edge_diff_rev <- edge_diff %>%
  rename(
    node1_tmp = node1,
    node2_tmp = node2
  ) %>%
  transmute(
    node1 = node2_tmp,
    node2 = node1_tmp,
    mean_diff = mean_diff,
    p = p,
    CI_lower = CI_lower,
    CI_upper = CI_upper,
    padj = padj
  )

edge_diff_bidir <- bind_rows(edge_diff, edge_diff_rev)

organs <- sort(unique(c(edge_diff_bidir$node1, edge_diff_bidir$node2)))

full_pairs <- expand.grid(
  node1 = organs,
  node2 = organs,
  stringsAsFactors = FALSE
) %>%
  filter(node1 != node2)

edge_full <- full_pairs %>%
  left_join(
    edge_diff_bidir %>%
      select(node1, node2, mean_diff, p, padj),
    by = c("node1", "node2")
  ) %>%
  mutate(
    diff = ifelse(is.na(mean_diff), 0, mean_diff),
    p    = ifelse(is.na(p), 1, p),
    padj = ifelse(is.na(padj), 1, padj)
  )

make_matrix <- function(edge_full, value_col, organs) {
  edge_full %>%
    mutate(
      node1 = factor(node1, levels = organs),
      node2 = factor(node2, levels = organs)
    ) %>%
    select(node1, node2, {{ value_col }}) %>%
    pivot_wider(
      names_from = node2,
      values_from = {{ value_col }}
    ) %>%
    arrange(node1) %>%
    column_to_rownames("node1") %>%
    as.matrix()
}

diff_mat <- make_matrix(edge_full, diff, organs)
p_mat    <- make_matrix(edge_full, p, organs)
padj_mat <- make_matrix(edge_full, padj, organs)

diff_mat <- (diff_mat + t(diff_mat)) / 2
p_mat    <- (p_mat + t(p_mat)) / 2
padj_mat <- (padj_mat + t(padj_mat)) / 2

diag(diff_mat) <- 0
diag(p_mat) <- 1
diag(padj_mat) <- 1

## 9. Export results
write.csv(diff_mat, file.path(outdir, "diff_edges_mat.csv"))
write.csv(p_mat, file.path(outdir, "diff_edges_p.csv"))
write.csv(padj_mat, file.path(outdir, "diff_edges_padj.csv"))

## 10. Plotting
diff_mat_scaled <- diff_mat / 0.2

corrplot(diff_mat_scaled,
  p.mat = p_mat,
  method = "circle",
  order = "AOE",
  type = "lower",
  diag = FALSE,
  insig = "label_sig", sig.level = c(0.01, 0.05, 0.10), pch.cex = 0.5,
  col = COL2("RdBu", 8)
)

col2 <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(6)

corrplot(as.matrix(diff_mat_scaled),
  p.mat = as.matrix(p_mat),
  order = "FPC",
  type = "lower",
  diag = FALSE,
  insig = "label_sig", sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.5,
  col = col2
)
}
