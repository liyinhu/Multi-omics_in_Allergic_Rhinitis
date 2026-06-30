# Usage:
#   Rscript "01.DIABOLO-model.R"
#   Rscript "01.DIABOLO-model.R" --proteome proteome.csv --metabolome metabolome.csv --lipidome lipidome.csv --group group.csv --outdir results

library(mixOmics)

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

read_matrix_csv <- function(path, label) {
  if (!file.exists(path)) {
    stop(sprintf("%s not found: %s", label, path))
  }
  read.csv(path, header = TRUE, row.names = 1, check.names = FALSE)
}

read_group_table <- function(path) {
  if (!file.exists(path)) {
    stop(sprintf("Group file not found: %s", path))
  }
  read.csv(path, header = TRUE, row.names = 1, check.names = FALSE)
}

align_to_common_samples <- function(proteome, metabolome, lipidome, group_df) {
  sample_sets <- list(
    rownames(proteome),
    rownames(metabolome),
    rownames(lipidome),
    rownames(group_df)
  )
  common_samples <- Reduce(intersect, sample_sets)
  if (length(common_samples) == 0) {
    stop("No overlapping sample IDs were found across the omics matrices and group table.")
  }
  proteome <- proteome[common_samples, , drop = FALSE]
  metabolome <- metabolome[common_samples, , drop = FALSE]
  lipidome <- lipidome[common_samples, , drop = FALSE]
  group_df <- group_df[common_samples, , drop = FALSE]
  list(proteome = proteome, metabolome = metabolome, lipidome = lipidome, group_df = group_df)
}

args <- commandArgs(trailingOnly = TRUE)
proteome_path <- parse_cli_value(args, c("--proteome"))
metabolome_path <- parse_cli_value(args, c("--metabolome"))
lipidome_path <- parse_cli_value(args, c("--lipidome"))
group_path <- parse_cli_value(args, c("--group"))
outdir <- parse_cli_value(args, c("--outdir"))

if (is.null(outdir) || !nzchar(outdir)) {
  outdir <- "."
}
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

if (is.null(proteome_path)) {
  proteome_path <- file.choose()
}
if (is.null(metabolome_path)) {
  metabolome_path <- file.choose()
}
if (is.null(lipidome_path)) {
  lipidome_path <- file.choose()
}
if (is.null(group_path)) {
  group_path <- file.choose()
}

# Import proteome, metabolome and lipidome
proteome <- read_matrix_csv(proteome_path, "Proteome matrix")
metabolome <- read_matrix_csv(metabolome_path, "Metabolome matrix")
lipidome <- read_matrix_csv(lipidome_path, "Lipidome matrix")
Y <- read_group_table(group_path)

if (nrow(proteome) == 0 || ncol(proteome) == 0) {
  stop("Proteome matrix is empty.")
}
if (nrow(metabolome) == 0 || ncol(metabolome) == 0) {
  stop("Metabolome matrix is empty.")
}
if (nrow(lipidome) == 0 || ncol(lipidome) == 0) {
  stop("Lipidome matrix is empty.")
}
if (nrow(Y) == 0 || ncol(Y) == 0) {
  stop("Group table is empty.")
}

if (is.null(rownames(proteome)) || is.null(rownames(metabolome)) || is.null(rownames(lipidome)) || is.null(rownames(Y))) {
  stop("All inputs must have sample IDs in row names.")
}

if (!("Group" %in% colnames(Y))) {
  stop("Group table must contain a Group column.")
}

aligned <- align_to_common_samples(proteome, metabolome, lipidome, Y)
proteome <- aligned$proteome
metabolome <- aligned$metabolome
lipidome <- aligned$lipidome
Y <- aligned$group_df

if (length(unique(na.omit(Y$Group))) < 2) {
  stop("At least two groups are required in the group table after alignment.")
}

data = list(Proteome = proteome,
            Metabolome = metabolome,
            Lipidome = lipidome)

# Pairwise PLS Comparisons
# Generate three pairwise PLS models
list.keepX = c(25, 25)
list.keepY = c(25, 25)

pls1 <- spls(data[["Metabolome"]], data[["Lipidome"]],
             keepX = list.keepX, keepY = list.keepY)
pls2 <- spls(data[["Proteome"]], data[["Metabolome"]],
             keepX = list.keepX, keepY = list.keepY)
pls3 <- spls(data[["Proteome"]], data[["Lipidome"]],
             keepX = list.keepX, keepY = list.keepY)


# Plot features of PLS
plotVar(pls1, cutoff = 0.5, title = "(a) Metabolome vs Lipidome",
        legend = c("Metabolome", "Lipidome"),
        var.names = FALSE, style = 'graphics',
        pch = c(16, 17), cex = c(2,2),
        col = c('#D86C76', '#14976F'))


plotVar(pls2, cutoff = 0.5, title = "(b) Proteome vs Metabolome",
        legend = c("Proteome", "Metabolome"),
        var.names = FALSE, style = 'graphics',
        pch = c(16, 17), cex = c(2,2),
        col = c('#58ABDB', '#D86C76'))


plotVar(pls3, cutoff = 0.5, title = "(c) Proteome vs Lipidome",
        legend = c("Proteome", "Lipidome"),
        var.names = FALSE, style = 'graphics',
        pch = c(16, 17), cex = c(2,2),
        col = c('#58ABDB', '#14976F'))

# Initial DIABLO Model
# For square matrix filled with 0.1s
design = matrix(0.1, ncol = length(data), nrow = length(data),
                dimnames = list(names(data), names(data)))
diag(design) = 0 # set diagonal to 0s

# Form basic DIABLO model
basic.diablo.model = block.splsda(X = data, Y = Y, ncomp = 5, design = design)


# Tuning the number of components
# Run component number tuning with repeated CV
perf.diablo = perf(basic.diablo.model, validation = 'Mfold',
                   folds = 10, nrepeat = 10)

plot(perf.diablo) # plot output of tuning

# Set the optimal ncomp value
ncomp = perf.diablo$choice.ncomp$WeightedVote["Overall.BER", "centroids.dist"]
# Show the optimal choice for ncomp for each dist metric
perf.diablo$choice.ncomp$WeightedVote

# Tuning the number of features
# Set grid of values for each component to test
test.keepX = list (mRNA = c(5:9, seq(10, 18, 2), seq(20,30,5)),
                   miRNA = c(5:9, seq(10, 18, 2), seq(20,30,5)),
                   proteomics = c(5:9, seq(10, 18, 2), seq(20,30,5)))

if (!all(names(test.keepX) %in% names(data))) {
  warning(sprintf(
    "test.keepX block names do not exactly match the data blocks: %s",
    paste(setdiff(names(test.keepX), names(data)), collapse = ", ")
  ))
}

# Run the feature selection tuning
tune.TCGA = tune.block.splsda(X = data, Y = Y, ncomp = ncomp,
                              test.keepX = test.keepX, design = design,
                              validation = 'Mfold', folds = 10, nrepeat = 1,
                              dist = "centroids.dist")

# Final DIABLO Model
# Set the optimised DIABLO model
final.diablo.model = block.splsda(X = data, Y = Y, ncomp = ncomp,
                          keepX = list.keepX, design = design)

final.diablo.model$design #Design matrix for the final model


# Plots
plot_path <- file.path(outdir, "DIABLO_plot.pdf")
circos_path <- file.path(outdir, "DIABLO_circos.pdf")
pdf(plot_path, width = 8, height = 6)
plotDiablo(final.diablo.model, ncomp = 1)
dev.off()

pdf(circos_path, width = 8, height = 6)
circosPlot(final.diablo.model, cutoff = 0.7, line = TRUE,
           color.blocks= c('darkorchid', 'brown1', 'lightgreen'),
           color.cor = c("chocolate3","grey20"), size.labels = 1.5)
dev.off()

# Export multi-omics relations
net <- network(final.diablo.model, blocks = c(1,2,3),
        color.node = c('darkorchid', 'brown1', 'lightgreen'), cutoff = 0.4)
edges <- net$gR$edges
write.csv(edges, file.path(outdir, "DIABLO_network_edges.csv"))
