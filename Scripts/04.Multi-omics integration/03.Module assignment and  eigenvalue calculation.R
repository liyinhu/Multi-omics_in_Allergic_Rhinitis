# Usage:
#   Rscript "03.Module assignment and  eigenvalue calculation.R"
#   Rscript "03.Module assignment and  eigenvalue calculation.R" --diablo-edges edges_from_DIABOLO.csv --mmvec-edges edges_from_mmVec.csv --signature multi_omics_signatures.csv --feature-type feature_type.csv --outdir results

library(MetaNet)
library(dplyr)
library(tibble)
library(WGCNA)
library(tidyr)

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

normalize_edge_table <- function(edge_df, source_label) {
  edge_df <- as.data.frame(edge_df, check.names = FALSE)
  if (all(c("Source", "Target") %in% names(edge_df))) {
    return(edge_df)
  }
  if (all(c("node1", "node2") %in% names(edge_df))) {
    names(edge_df)[match(c("node1", "node2"), names(edge_df))] <- c("Source", "Target")
    message(sprintf("Mapped %s edge columns node1/node2 to Source/Target.", source_label))
    return(edge_df)
  }
  if (all(c("from", "to") %in% names(edge_df))) {
    names(edge_df)[match(c("from", "to"), names(edge_df))] <- c("Source", "Target")
    message(sprintf("Mapped %s edge columns from/to to Source/Target.", source_label))
    return(edge_df)
  }
  if (ncol(edge_df) >= 2) {
    names(edge_df)[1:2] <- c("Source", "Target")
    warning(sprintf(
      "%s does not contain Source/Target or node1/node2 columns; the first two columns were treated as Source/Target.",
      source_label
    ))
    return(edge_df)
  }
  stop(sprintf("%s must contain at least two columns for Source and Target.", source_label))
}

read_signature_matrix <- function(path) {
  sig <- read.csv(path, header = TRUE, row.names = 1, check.names = FALSE)
  if (nrow(sig) == 0 || ncol(sig) == 0) {
    stop("signature file is empty.")
  }
  sig_numeric <- suppressWarnings(as.data.frame(lapply(sig, function(x) as.numeric(as.character(x))), check.names = FALSE))
  bad_cols <- names(sig)[vapply(seq_along(sig), function(i) any(is.na(sig_numeric[[i]]) & !is.na(sig[[i]])), logical(1))]
  if (length(bad_cols) > 0) {
    stop(sprintf("signature file contains non-numeric values in columns: %s", paste(bad_cols, collapse = ", ")))
  }
  sig_numeric
}

read_feature_type <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("tsv", "txt")) {
    ft <- read.delim(path, header = TRUE, check.names = FALSE)
  } else {
    ft <- read.csv(path, header = TRUE, check.names = FALSE)
  }
  required_cols <- c("Feature", "Type")
  missing_cols <- setdiff(required_cols, names(ft))
  if (length(missing_cols) > 0) {
    stop(sprintf("feature_type file is missing required columns: %s", paste(missing_cols, collapse = ", ")))
  }
  ft <- ft[, required_cols, drop = FALSE]
  ft$Feature <- as.character(ft$Feature)
  ft$Type <- as.character(ft$Type)
  ft
}

args <- commandArgs(trailingOnly = TRUE)
diablo_edges_path <- parse_cli_value(args, c("--diablo-edges"))
mmvec_edges_path <- parse_cli_value(args, c("--mmvec-edges"))
signature_path <- parse_cli_value(args, c("--signature"))
feature_type_path <- parse_cli_value(args, c("--feature-type"))
outdir <- parse_cli_value(args, c("--outdir"))

if (is.null(diablo_edges_path)) {
  diablo_edges_path <- "edges_from_DIABOLO.csv"
}
if (is.null(mmvec_edges_path)) {
  mmvec_edges_path <- "edges_from_mmVec.csv"
}
if (is.null(signature_path)) {
  signature_path <- "multi_omics_signatures.csv"
}
if (is.null(feature_type_path)) {
  feature_type_path <- "feature_type.csv"
}
if (is.null(outdir) || !nzchar(outdir)) {
  outdir <- "."
}

if (!file.exists(diablo_edges_path)) {
  stop(sprintf("DIABLO edge file not found: %s", diablo_edges_path))
}
if (!file.exists(mmvec_edges_path)) {
  stop(sprintf("mmVec edge file not found: %s", mmvec_edges_path))
}
if (!file.exists(signature_path)) {
  stop(sprintf("signature file not found: %s", signature_path))
}
if (!file.exists(feature_type_path)) {
  stop(sprintf("feature_type file not found: %s", feature_type_path))
}

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Import edges from DIABOLO and mmVec models
edge_list1 <- normalize_edge_table(read.csv(diablo_edges_path, header = TRUE, check.names = FALSE), "DIABOLO")
edge_list2 <- normalize_edge_table(read.csv(mmvec_edges_path, header = TRUE, check.names = FALSE), "mmVec")

if (!all(c("Source", "Target") %in% names(edge_list1)) || !all(c("Source", "Target") %in% names(edge_list2))) {
  stop("Both edge files must provide Source and Target columns after normalization.")
}

edge_list <- bind_rows(edge_list1, edge_list2)

if (nrow(edge_list) == 0) {
  stop("Combined edge list is empty.")
}

if (any(is.na(edge_list$Source)) || any(is.na(edge_list$Target))) {
  stop("Edge list contains missing Source or Target values.")
}

# Construct multi-omics network
dnet <- c_net_from_edgelist(edge_list, direct = F)

# Module detection
co_net_modu <- module_detect(dnet, method = "cluster_fast_greedy")
node_modu <- get_v(co_net_modu)[, c("name", "module")]
write.csv(node_modu, file.path(outdir, "Node_module.csv"), row.names = FALSE)

# Plot
pdf(file.path(outdir, "Module_tree.pdf"), width = 8, height = 6)
plot_module_tree(co_net_modu, label.size = 0.6)
dev.off()

## Covert to long format: rows and columns represent features and samples, respectively
expr_mat <- read_signature_matrix(signature_path)
feature_type <- read_feature_type(feature_type_path)

all_features <- unique(c(rownames(expr_mat), edge_list$Source, edge_list$Target))
missing_in_type <- setdiff(all_features, feature_type$Feature)
if (length(missing_in_type) > 0) {
  warning(sprintf(
    "%d features appearing in edges/signatures are not present in feature_type; they will be kept with Type = NA.",
    length(missing_in_type)
  ))
}

expr_long <- as.data.frame(t(expr_mat), check.names = FALSE)
expr_long$Feature <- rownames(expr_long)

# Omics integration
expr_annot <- expr_long %>%
  left_join(feature_type, by = "Feature")

# Z-score normalization by features
expr_scaled_long <- expr_annot %>%
  group_by(Type) %>%
  mutate(across(where(is.numeric), scale)) %>%
  ungroup()

# Convert to wide format
expr_scaled <- expr_scaled_long %>%
  select(-Type) %>%
  column_to_rownames("Feature") %>%
  as.data.frame()

expr_scaled <- t(expr_scaled)

node_modu$module <- as.character(node_modu$module)

# Ensure feature consistency
expr_scaled <- expr_scaled[intersect(rownames(expr_scaled), node_modu$name), ]
node_modu <- node_modu %>% filter(name %in% rownames(expr_scaled))

if (nrow(expr_scaled) == 0 || nrow(node_modu) == 0) {
  stop("No overlapping features were found between the signature matrix and the detected modules.")
}

# Calculate the eigengene by modules
modules_list <- split(node_modu$name, node_modu$module)

datExpr <- t(expr_scaled)

# Calcualte module eigengenes
ME_list <- lapply(names(modules_list), function(mod) {
  subset_features <- modules_list[[mod]]
  if (length(subset_features) >= 2) {
    dat_mod <- datExpr[, subset_features, drop = FALSE]

    # Process NA
    dat_mod <- as.matrix(dat_mod)
    storage.mode(dat_mod) <- "numeric"
    dat_mod <- dat_mod[, colSums(is.na(dat_mod)) == 0, drop = FALSE]
    
    # moduleEigengenes
    if (ncol(dat_mod) >= 2) {
      ME <- moduleEigengenes(dat_mod, colors = rep(mod, ncol(dat_mod)))$eigengenes
      colnames(ME) <- mod
      return(ME)
    }
  }
  return(NULL)
})

ME_list <- Filter(Negate(is.null), ME_list)
if (length(ME_list) == 0) {
  stop("No modules satisfied the minimum feature requirement for moduleEigengenes.")
}

# Export the results
ME_matrix <- do.call(cbind, ME_list)
write.csv(ME_matrix, file.path(outdir, "ME_matrix.csv"))
