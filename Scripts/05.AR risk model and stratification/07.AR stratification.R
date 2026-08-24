# Usage:
#   Rscript "07.AR stratification.R"
#   Rscript "07.AR stratification.R" --input target_matrix.csv --outdir results

library(cluster)
library(factoextra)
library(dplyr)
library(tidyr)
library(FSA)
library(ggplot2)
library(umap)
library(vegan)
library(RColorBrewer)

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

if (is.null(outdir) || !nzchar(outdir)) {
  outdir <- "."
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

# Import target multi-omics data
data <- read.csv(input_file, header = TRUE, row.names = 1, check.names = FALSE)

if (nrow(data) == 0 || ncol(data) == 0) {
  stop("Input matrix is empty.")
}
if (is.null(rownames(data)) || any(!nzchar(rownames(data)))) {
  stop("Sample IDs are missing from the row names.")
}

data_numeric <- suppressWarnings(data.frame(lapply(data, function(x) as.numeric(as.character(x))), check.names = FALSE))
non_numeric_cols <- names(data)[vapply(seq_along(data), function(i) any(is.na(data_numeric[[i]]) & !is.na(data[[i]])), logical(1))]
if (length(non_numeric_cols) > 0) {
  stop(sprintf("The following feature columns cannot be converted to numeric: %s", paste(non_numeric_cols, collapse = ", ")))
}
data <- data_numeric

# === Step 1: Determine the optimal number of clusters ===
max_k <- 10
wss <- numeric(max_k)
wss[1] <- sum(scale(data, center = TRUE, scale = FALSE)^2)

set.seed(666)
for (k in 2:max_k) {
  km <- kmeans(data, centers = k, nstart = 5)
  wss[k] <- km$tot.withinss
}

elbow_plot <- ggplot(data.frame(k = 1:max_k, wss = wss), aes(x = k, y = wss)) +
  geom_line(color = "#58ABDB") +
  geom_point(color = "#58ABDB") +
  geom_vline(xintercept = which.min(diff(diff(wss))), color = "red", linetype = 2) +
  labs(
    x = "Number of Clusters (K)",
    y = "Within-cluster sum of squares",
    title = "Elbow Curve"
  ) +
  theme_bw()

print(elbow_plot)
ggsave(file.path(outdir, "Elbow_Curve.pdf"), elbow_plot, width = 7, height = 5)

# === Step 2: Clustering (Optimal k can be adjusted according to the above figure)===
best_k <- 5
set.seed(123)
km_res <- kmeans(data, centers = best_k, nstart = 10)

data_clustered <- data.frame(Sample = rownames(data),
                             Cluster = factor(km_res$cluster),
                             data)

# === Step 3: PCoA visualization ===

dist_mat <- dist(data)

pcoa_res <- cmdscale(dist_mat, k = 2, eig = TRUE)

pcoa_df <- data.frame(
  Sample = rownames(data),
  PCoA1 = pcoa_res$points[, 1],
  PCoA2 = pcoa_res$points[, 2],
  Cluster = factor(km_res$cluster)
)

system_colors <- structure(
    RColorBrewer::brewer.pal(length(unique(km_res$cluster)), "Set2"),
    names = unique(km_res$cluster)
)

eig <- pcoa_res$eig

pcoa_plot <- ggplot(pcoa_df) +
  geom_point(aes(x = PCoA1, y = PCoA2, color = Cluster)) +
  stat_ellipse(aes(x = PCoA1, y = PCoA2, group = Cluster, color = Cluster), level = 0.95) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  labs(
    x = paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits = 4), "%)", sep = ""),
    y = paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits = 4), "%)", sep = "")
  ) +
  scale_color_manual(values = system_colors)

print(pcoa_plot)
ggsave(file.path(outdir, "PCoA.pdf"), pcoa_plot, width = 7, height = 5)

# === Step 4: Remove clusters with ≥ 5 samples ===
min_cluster_size <- 5
cluster_size <- table(data_clustered$Cluster)
valid_clusters <- names(cluster_size[cluster_size >= min_cluster_size])

if (length(valid_clusters) == 0) {
  stop("No clusters meet the minimum cluster size requirement.")
}

filtered_df <- data_clustered %>% filter(Cluster %in% valid_clusters)

if (length(unique(filtered_df$Cluster)) < 2) {
  stop("Fewer than two clusters remain after applying the minimum cluster size filter.")
}

long_df <- filtered_df %>%
  pivot_longer(cols = -c(Sample, Cluster), names_to = "Feature", values_to = "Value")

if (nrow(long_df) == 0) {
  stop("No data remain after filtering clusters by minimum cluster size.")
}

# === Step 5: Kruskal-Wallis + DunnTest Analysis ===
diff_res <- long_df %>%
  group_by(Feature) %>%
  summarise(p = kruskal.test(Value ~ Cluster)$p.value, .groups = "drop") %>%
  mutate(p_adj = p.adjust(p, method = "BH")) %>%
  arrange(p_adj)

sig_feats <- diff_res %>% filter(p_adj < 0.05) %>% pull(Feature)

if (length(sig_feats) == 0) {
  message("No significant features were found by Kruskal-Wallis after BH adjustment. Writing empty downstream results.")
  dunn_all <- data.frame()
  rep_features <- data.frame()
  top_cluster_feats <- data.frame()
} else {
  dunn_results <- list()
  for (f in sig_feats) {
    subset_data <- long_df %>% filter(Feature == f)
    dunn <- dunnTest(Value ~ Cluster, data = subset_data, method = "bh")
    dunn_df <- dunn$res
    dunn_df$Feature <- f
    dunn_results[[f]] <- dunn_df
  }
  if (length(dunn_results) == 0) {
    dunn_all <- data.frame()
  } else {
    dunn_all <- do.call(rbind, dunn_results)
  }

  # === Step 6: Isolate representative features for the clusters ===
  cluster_means <- long_df %>%
    filter(Feature %in% sig_feats) %>%
    group_by(Feature, Cluster) %>%
    summarise(mean_val = mean(Value), .groups = "drop")

  top_cluster_feats <- cluster_means %>%
    group_by(Feature) %>%
    filter(mean_val == max(mean_val)) %>%
    ungroup()

  rep_features <- top_cluster_feats %>%
    group_by(Cluster) %>%
    summarise(Representative_Features = paste(Feature, collapse = ", "), .groups = "drop")
}

print(rep_features)

# === Step 7: Visualize the representative features  ===
if (exists("top_cluster_feats") && nrow(top_cluster_feats) > 0) {
  for (f in unique(top_cluster_feats$Feature)) {
    p <- ggplot(long_df %>% filter(Feature == f),
                aes(x = Cluster, y = Value, fill = Cluster)) +
      geom_violin(trim = FALSE, alpha = 0.6) +
      geom_boxplot(width = 0.1, outlier.shape = NA, color = "black") +
      theme_bw() +
      scale_fill_manual(values = system_colors) +
      labs(title = paste("Violin Plot of", f)) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
    print(p)
    ggsave(file.path(outdir, paste0("Violin_", make.names(f), ".pdf")), p, width = 7, height = 5)
  }
}

# === Step 8: Export results ===
write.csv(diff_res, file.path(outdir, "Kruskal_Wallis_pvalues.csv"), row.names = FALSE)
write.csv(dunn_all, file.path(outdir, "DunnTest_pairwise_results.csv"), row.names = FALSE)
write.csv(rep_features, file.path(outdir, "Cluster_Representative_Features.csv"), row.names = FALSE)
write.csv(data_clustered, file.path(outdir, "Sample_clustered.csv"), row.names = FALSE)
