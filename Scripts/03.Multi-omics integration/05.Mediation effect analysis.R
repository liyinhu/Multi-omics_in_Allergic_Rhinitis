# Usage:
#   Rscript "05.Mediation effect analysis.R"
#   Rscript "05.Mediation effect analysis.R" --subpath Subpath.csv --signature multi_omics_signatures.csv --outdir results

library(lme4)
library(mediation)
library(dplyr)

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
subpath_file <- parse_cli_value(args, c("--subpath"))
signature_file <- parse_cli_value(args, c("--signature"))
outdir <- parse_cli_value(args, c("--outdir"))

if (is.null(subpath_file)) {
  subpath_file <- "Subpath.csv"
}
if (is.null(signature_file)) {
  signature_file <- "multi_omics_signatures.csv"
}
if (is.null(outdir) || !nzchar(outdir)) {
  outdir <- "."
}

if (!file.exists(subpath_file)) {
  stop(sprintf("Subpath file not found: %s", subpath_file))
}
if (!file.exists(signature_file)) {
  stop(sprintf("Signature file not found: %s", signature_file))
}

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Import isolated subpath
triplets_df <- read.csv(subpath_file, header = TRUE, check.names = FALSE)
all_data <- read.csv(signature_file, header = TRUE, row.names = 1, check.names = FALSE)

required_triplet_cols <- c("Source", "Middle", "Target")
missing_triplet_cols <- setdiff(required_triplet_cols, names(triplets_df))
if (length(missing_triplet_cols) > 0) {
    stop(sprintf("Subpath file is missing required columns: %s", paste(missing_triplet_cols, collapse = ", ")))
}

if (!"Group" %in% names(all_data)) {
    stop("Signature file must contain a Group column.")
}

# Ensure the requested features exist
required_features <- unique(c(triplets_df$Source, triplets_df$Middle, triplets_df$Target))
missing_features <- setdiff(required_features, names(all_data))
if (length(missing_features) > 0) {
    warning(sprintf(
        "The following features are missing from the signature file and corresponding paths will be skipped: %s",
        paste(missing_features, collapse = ", ")
    ))
}

results_list <- list()

# 双向中介效应分析
for (i in seq_len(nrow(triplets_df))) {
    microbe <- triplets_df$Source[i]
    metabolite <- triplets_df$Middle[i]
    protein <- triplets_df$Target[i]

    if (any(!c(microbe, metabolite, protein) %in% names(all_data))) {
        warning(sprintf("Skipping path %s | %s | %s because one or more features are missing from the signature file.", microbe, metabolite, protein))
        next
    }

    dat <- data.frame(
        microbe = all_data[[microbe]],
        metabolite = all_data[[metabolite]],
        protein = all_data[[protein]],
        Group = all_data$Group,
        check.names = FALSE
    )

    numeric_cols <- c("microbe", "metabolite", "protein")
    non_numeric <- numeric_cols[!vapply(dat[numeric_cols], is.numeric, logical(1))]
    if (length(non_numeric) > 0) {
        warning(sprintf(
            "Skipping path %s | %s | %s because variables are not numeric: %s",
            microbe, metabolite, protein, paste(non_numeric, collapse = ", ")
        ))
        next
    }

    if (sum(complete.cases(dat[, numeric_cols])) < 3) {
        warning(sprintf(
            "Skipping path %s | %s | %s because of too many missing values.",
            microbe, metabolite, protein
        ))
        next
    }

    # Model A：microbe → metabolite → protein
    model.m1 <- tryCatch(lm(metabolite ~ microbe, dat), error = function(e) e)
    model.y1 <- tryCatch(lm(protein ~ microbe + metabolite, data = dat), error = function(e) e)
    if (inherits(model.m1, "error") || inherits(model.y1, "error")) {
        warning(sprintf(
            "Skipping path %s | %s | %s because model fitting failed in direction A.",
            microbe, metabolite, protein
        ))
        next
    }

    med1 <- tryCatch(
        mediate(model.m1, model.y1, treat = "microbe", mediator = "metabolite", sims = 500),
        error = function(e) e
    )
    if (inherits(med1, "error")) {
        warning(sprintf(
            "Skipping path %s | %s | %s because mediation failed in direction A: %s",
            microbe, metabolite, protein, med1$message
        ))
        next
    }

    # Model B：protein → metabolite → microbe
    model.m2 <- tryCatch(lm(metabolite ~ protein, dat), error = function(e) e)
    model.y2 <- tryCatch(lm(microbe ~ protein + metabolite, dat), error = function(e) e)
    if (inherits(model.m2, "error") || inherits(model.y2, "error")) {
        warning(sprintf(
            "Skipping path %s | %s | %s because model fitting failed in direction B.",
            microbe, metabolite, protein
        ))
        next
    }

    med2 <- tryCatch(
        mediate(model.m2, model.y2, treat = "protein", mediator = "metabolite", sims = 500),
        error = function(e) e
    )
    if (inherits(med2, "error")) {
        warning(sprintf(
            "Skipping path %s | %s | %s because mediation failed in direction B: %s",
            microbe, metabolite, protein, med2$message
        ))
        next
    }

    # Results
    extract_summary <- function(med, direction) {
        data.frame(
            Direction = direction,
            ACME = med$d0,
            ACME_CI_Lower = med$d0.ci[1],
            ACME_CI_Upper = med$d0.ci[2],
            ACME_p = med$d0.p,
            
            ADE = med$z0,
            ADE_CI_Lower = med$z0.ci[1],
            ADE_CI_Upper = med$z0.ci[2],
            ADE_p = med$z0.p,
            
            Total_Effect = med$tau.coef,
            Total_CI_Lower = med$tau.ci[1],
            Total_CI_Upper = med$tau.ci[2],
            Total_p = med$tau.p,
            
            Prop_Mediated = med$n0,
            Prop_p = med$n0.p
        )
    }

    res_A <- extract_summary(med1, "Microbe鈫扢etabolite鈫扨rotein")
    res_B <- extract_summary(med2, "Protein鈫扢etabolite鈫扢icrobe")

    result_row <- bind_rows(res_A, res_B) %>%
        mutate(
            Microbe = microbe,
            Metabolite = metabolite,
            Protein = protein,
            Triplet = paste(microbe, metabolite, protein, sep = " | ")
        )

    results_list[[length(results_list) + 1]] <- result_row
}

if (length(results_list) == 0) {
    warning("No mediation results were generated; writing an empty Mediation_results.csv.")
    final_results <- data.frame()
} else {
    final_results <- bind_rows(results_list)
}

write.csv(final_results, file.path(outdir, "Mediation_results.csv"), row.names = FALSE)
