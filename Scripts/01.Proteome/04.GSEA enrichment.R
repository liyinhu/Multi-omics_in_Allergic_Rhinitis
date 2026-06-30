# Usage:
#   Rscript "04.GSEA enrichment.R"
#   Rscript "04.GSEA enrichment.R" --input Proteome_summary.csv --outdir results --protein-column ProteinID --estimate-column AR_Estimate

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(forcats)
library(dplyr)
library(ggplot2)

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

args <- commandArgs(trailingOnly = TRUE)
input_path <- parse_cli_value(args, c("--input", "-i"))
outdir <- parse_cli_value(args, c("--outdir", "-o"))
protein_column <- parse_cli_value(args, c("--protein-column", "-p"))
estimate_column <- parse_cli_value(args, c("--estimate-column", "-e"))

if (is.null(outdir) || !nzchar(outdir)) {
  outdir <- "."
}
if (is.null(protein_column) || !nzchar(protein_column)) {
  protein_column <- "ProteinID"
}
if (is.null(estimate_column) || !nzchar(estimate_column)) {
  estimate_column <- "AR_Estimate"
}

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Import differential protein analysis results
if (is.null(input_path) || !nzchar(input_path)) {
  input_path <- file.choose()
}

if (!file.exists(input_path)) {
  stop(sprintf("Input file not found: %s", input_path))
}

data <- read.csv(input_path, header = TRUE, check.names = FALSE)

if (!(protein_column %in% names(data))) {
  if (identical(protein_column, "ProteinID") && ncol(data) >= 1) {
    warning(sprintf(
      "ProteinID column not found; using the first column as protein IDs: %s",
      names(data)[1]
    ))
    protein_column <- names(data)[1]
  } else {
    stop(sprintf("Protein ID column not found: %s", protein_column))
  }
}

if (!(estimate_column %in% names(data))) {
  stop(sprintf("Estimate column not found: %s", estimate_column))
}

estimate_values <- suppressWarnings(as.numeric(as.character(data[[estimate_column]])))
if (any(is.na(estimate_values) & !is.na(data[[estimate_column]]))) {
  stop(sprintf("Estimate column must be numeric or convertible to numeric: %s", estimate_column))
}

# Isolate protein ID and estimate size
df <- data.frame(
  ProteinID = as.character(data[[protein_column]]),
  AR_Estimate = estimate_values,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
df$SYMBOL <- df$ProteinID

df <- df[!is.na(df$ProteinID) & nzchar(df$ProteinID) & !is.na(df$AR_Estimate), , drop = FALSE]
if (nrow(df) == 0) {
  stop("No valid protein IDs and AR_Estimate values remain after removing missing values.")
}

# Protein ID transformation
entrez <- bitr(
  df$ProteinID,
  fromType = "SYMBOL",
  toType = c("ENTREZID"),
  OrgDb = "org.Hs.eg.db"
)
df <- merge(df, entrez, by.y = "SYMBOL")

# Export gene list ordering input
geneList <- df$AR_Estimate
names(geneList) <- df$ENTREZID
geneList <- sort(geneList, decreasing = TRUE)
geneList <- geneList[!is.na(geneList)]

if (length(geneList) == 0) {
  stop("geneList is empty after ID conversion and NA removal.")
}

# GSEA KEGG enrichment analysis
kegg <- gseKEGG(
  geneList,
  organism = "hsa",
  pvalueCutoff = 0.2,
  pAdjustMethod = "BH",
  minGSSize = 5,
  maxGSSize = 200
)
kegg <- append_kegg_category(kegg)

# Convert to readable mode
kegg <- setReadable(
  kegg,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID"
)

# Export enrichment results
kegg_result <- kegg@result
write.csv(kegg_result, file.path(outdir, "AR_gsea_KEGG.csv"), row.names = FALSE)

# Plotting
kegg_result2 <- kegg_result %>% mutate(Description = fct_reorder(Description, NES))

kegg_plot <- ggplot(data = kegg_result2, aes(x = NES, y = Description)) +
  geom_point(aes(alpha = 0.7, size = -log10(pvalue), color = NES)) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(color = "white"),
    panel.grid.minor = element_line(color = "white"),
    legend.title = element_blank()
  ) +
  scale_color_gradient2(low = "#205D9B", high = "#D4423A", limits = c(-1.8, 1.8)) +
  labs(x = "Normalized Enrichment Score", y = "", title = "AR vs HC") +
  scale_size("-log10(P-value)", range = c(3, 7))

ggsave(
  filename = file.path(outdir, "AR_gsea_KEGG_dotplot.pdf"),
  plot = kegg_plot,
  width = 8,
  height = 6
)
