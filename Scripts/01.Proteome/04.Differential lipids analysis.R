# Usage:
#   Rscript "04.Differential lipids analysis.R"
#   Rscript "04.Differential lipids analysis.R" --input adjusted_lipidome.csv --outdir results --group-column Group

library(ggrepel)
library(ggplot2)
library(grid)

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
group_column <- parse_cli_value(args, c("--group-column", "-g"))

if (is.null(outdir) || !nzchar(outdir)) {
  outdir <- "."
}
if (is.null(group_column) || !nzchar(group_column)) {
  group_column <- "Group"
}

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Import adjusted lipidome
if (is.null(input_path) || !nzchar(input_path)) {
  input_path <- file.choose()
}

if (!file.exists(input_path)) {
  stop(sprintf("Input file not found: %s", input_path))
}

df <- read.csv(input_path, header = TRUE, check.names = FALSE)

if (!(group_column %in% names(df))) {
  stop(sprintf("Group column not found: %s", group_column))
}

if (nrow(df) == 0) {
  stop("Input file is empty.")
}

group_values <- df[[group_column]]
if (length(unique(na.omit(group_values))) < 2) {
  stop(sprintf("At least two groups are required in column %s.", group_column))
}

# Preserve the original compatibility assumption: the first two columns are non-lipid columns.
lipid_col_idx <- seq.int(3, ncol(df))
if (length(lipid_col_idx) == 0 || max(lipid_col_idx) > ncol(df)) {
  stop("At least one lipid feature column is required from the third column onward.")
}

lipid_df <- df[, lipid_col_idx, drop = FALSE]
if (ncol(lipid_df) == 0) {
  stop("No lipid feature columns were found from the third column onward.")
}

lipid_df_numeric <- suppressWarnings(
  as.data.frame(lapply(lipid_df, function(x) as.numeric(as.character(x))), check.names = FALSE)
)

bad_cols <- names(lipid_df)[vapply(
  seq_along(lipid_df),
  function(i) any(is.na(lipid_df_numeric[[i]]) & !is.na(lipid_df[[i]])),
  logical(1)
)]

if (length(bad_cols) > 0) {
  stop(sprintf(
    "Lipid feature columns must be numeric or convertible to numeric. Invalid columns: %s",
    paste(bad_cols, collapse = ", ")
  ))
}

df[, lipid_col_idx] <- lipid_df_numeric

Number_of_lipids <- ncol(lipid_df_numeric)

Test_summary <- data.frame(matrix(ncol = 4, nrow = Number_of_lipids))
colnames(Test_summary) <- c("AR_Estimate", "AR_SE", "AR_T_value", "AR_P_value")
rownames(Test_summary) <- names(df[c(3:(Number_of_lipids + 2))])

# Differential analysis between AR and HC using linear regression models
for (i in 1:Number_of_lipids) {
  tryCatch({
    formula_text <- sprintf("%s ~ %s", names(df)[i + 2], group_column)
    Test <- lm(stats::as.formula(formula_text), data = df)
    Test_summary[i, 1:4] <- as.vector(summary(Test)$coefficients[2, 1:4])
  }, error = function(e) {
    cat("ERROR :", conditionMessage(e), "\n")
  })
}

# BH adjustment for multiple testing
Test_summary$AR_FDR <- p.adjust(Test_summary$AR_P_value, method = "BH")

# Export statistical results
write.csv(Test_summary, file.path(outdir, "Lipidome_summary.csv"), row.names = TRUE)

# Volcano plot
k1 <- (Test_summary$AR_P_value < 0.05) & (Test_summary$AR_Estimate > 0)
k2 <- (Test_summary$AR_P_value < 0.05) & (Test_summary$AR_Estimate < 0)
Test_summary$change <- ifelse(k1, "UP", ifelse(k2, "DOWN", "NON"))
Test_summary$label <- ifelse(
  Test_summary$AR_P_value < 0.05 & (Test_summary$AR_Estimate < -0.575 | Test_summary$AR_Estimate > 0.575),
  as.character(rownames(Test_summary)),
  ""
)

volcano_plot <- ggplot(
  data = Test_summary,
  aes(x = AR_Estimate, y = -log10(AR_P_value), colour = change, fill = change)
) +
  geom_point(alpha = 0.6, aes(size = -log10(AR_P_value))) +
  geom_text_repel(
    aes(x = AR_Estimate, y = -log10(AR_P_value), label = label),
    size = 3,
    box.padding = unit(0.6, "lines"),
    point.padding = unit(0.7, "lines"),
    segment.color = "black",
    show.legend = FALSE
  ) +
  geom_hline(yintercept = -log10(0.05), color = "gray", size = 0.5) +
  theme_bw() +
  labs(
    x = "Normalized effect size",
    y = "-log10(P-value)",
    title = "Plasma lipids (AR vs HC)"
  ) +
  theme(
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 13),
    plot.title = element_text(hjust = 0.5, size = 15, face = "bold")
  ) +
  scale_color_manual(values = c("steelblue2", "gray30", "tomato2")) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(color = "black", fill = "transparent"),
    legend.title = element_blank(),
    legend.key = element_rect(fill = "transparent"),
    legend.background = element_rect(fill = "transparent")
  )

ggsave(
  filename = file.path(outdir, "Lipidome_volcano.pdf"),
  plot = volcano_plot,
  width = 8,
  height = 6
)
