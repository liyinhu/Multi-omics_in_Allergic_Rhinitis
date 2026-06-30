# Usage:
#   Rscript 04.Structural Equations Model.R
#   Rscript 04.Structural Equations Model.R --organ organ_matrix.csv --info info.tsv
#
# If no arguments are provided, the script falls back to interactive file selection.

library(lavaan)
library(semPlot)
library(psych)

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

common_samples <- intersect(rownames(organ), info$ID)
if (length(common_samples) == 0) {
  stop("No overlapping sample IDs were found between the organ matrix and sample info.")
}

organ <- organ[common_samples, , drop = FALSE]
info <- info[match(common_samples, info$ID), , drop = FALSE]

## 2. Z-score normalization
organ_z <- scale(organ)
organ_z <- as.data.frame(organ_z)

# 3. Group: NC = 0, AR = 1
info$Group_bin <- ifelse(info$Group == "2AR", 1, 0)

dat <- cbind(
  organ_z,
  Group = info$Group_bin
)

## 4. Construct SEM
model2 <- '
Lymphoidtissue ~ Bonemarrow
Lymphoidtissue ~ Gastrointestinal
Hepatobiliary ~ Lymphoidtissue
Hepatobiliary ~ Brain
Hepatobiliary ~ Endocrine
Brain ~ Gastrointestinal
Brain ~ Pancreas
Brain ~ Hepatobiliary
Brain ~ Endocrine
Gastrointestinal ~ Hepatobiliary
Gastrointestinal ~ Pancreas
Gastrointestinal ~ Endocrine
Bonemarrow ~ Brain
Bonemarrow ~ Gastrointestinal
Pancreas ~ Gastrointestinal
Endocrine ~ Lymphoidtissue
'

## 5. Model assessment
fit <- sem(model2, data = dat, estimator = "MLR", missing = "fiml", meanstructure = TRUE)
summary(fit, fit.measures = TRUE, standardized = TRUE)
