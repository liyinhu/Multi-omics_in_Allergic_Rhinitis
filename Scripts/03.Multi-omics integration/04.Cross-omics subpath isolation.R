# Usage:
#   Rscript "04.Cross-omics subpath isolation.R"
#   Rscript "04.Cross-omics subpath isolation.R" --edges multi-omics.network.edge --nodes multi-omics.network.nodes --outdir results

library(igraph)
library(dplyr)
library(tidyr)
library(purrr)

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

normalize_edge_table <- function(edge_df) {
  edge_df <- as.data.frame(edge_df, check.names = FALSE)
  lower_names <- tolower(names(edge_df))
  if (all(c("source", "target") %in% lower_names)) {
    names(edge_df)[match(c("source", "target"), lower_names)] <- c("Source", "Target")
    return(edge_df)
  }
  if (all(c("node1", "node2") %in% lower_names)) {
    names(edge_df)[match(c("node1", "node2"), lower_names)] <- c("Source", "Target")
    return(edge_df)
  }
  if (all(c("from", "to") %in% lower_names)) {
    names(edge_df)[match(c("from", "to"), lower_names)] <- c("Source", "Target")
    return(edge_df)
  }
  if (all(c("Source", "Target") %in% names(edge_df))) {
    return(edge_df)
  }
  stop("Edge file must contain Source/Target columns, or a compatible source/target, node1/node2, or from/to pair.")
}

normalize_node_table <- function(node_df) {
  node_df <- as.data.frame(node_df, check.names = FALSE)
  lower_names <- tolower(names(node_df))

  rename_map <- c(
    "id" = "ID",
    "module" = "module",
    "group" = "Group"
  )

  for (src in names(rename_map)) {
    if (src %in% lower_names) {
      names(node_df)[which(lower_names == src)] <- rename_map[[src]]
    }
  }

  node_df
}

args <- commandArgs(trailingOnly = TRUE)
edges_path <- parse_cli_value(args, c("--edges"))
nodes_path <- parse_cli_value(args, c("--nodes"))
outdir <- parse_cli_value(args, c("--outdir"))

if (is.null(edges_path)) {
  edges_path <- "multi-omics network.edge"
}
if (is.null(nodes_path)) {
  nodes_path <- "multi-omics network.nodes"
}
if (is.null(outdir) || !nzchar(outdir)) {
  outdir <- "."
}

if (!file.exists(edges_path)) {
  stop(sprintf("Edge file not found: %s", edges_path))
}
if (!file.exists(nodes_path)) {
  stop(sprintf("Node file not found: %s", nodes_path))
}

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Import data
edges <- normalize_edge_table(read.csv(edges_path, header = TRUE, check.names = FALSE))
nodes <- normalize_node_table(read.csv(nodes_path, header = TRUE, check.names = FALSE))

required_edge_cols <- c("Source", "Target")
missing_edge_cols <- setdiff(required_edge_cols, names(edges))
if (length(missing_edge_cols) > 0) {
  stop(sprintf("Edge file is missing required columns: %s", paste(missing_edge_cols, collapse = ", ")))
}

required_node_cols <- c("ID", "module", "Group")
missing_node_cols <- setdiff(required_node_cols, names(nodes))
if (length(missing_node_cols) > 0) {
  stop(sprintf("Node file is missing required columns: %s", paste(missing_node_cols, collapse = ", ")))
}

if (nrow(edges) == 0 || nrow(nodes) == 0) {
  stop("Edge or node table is empty.")
}

# Network construction
g <- graph_from_data_frame(d = edges, vertices = nodes, directed = FALSE)

# Module isolation
modules <- unique(nodes$module)

# Initialization
all_paths <- list()

for (mod in modules) {
    cat("Processing module:", mod, "\n")
    
    # Subpath isolation
    mod_nodes <- nodes %>% filter(module == mod)
    mod_edges <- edges %>% filter(Source %in% mod_nodes$ID & Target %in% mod_nodes$ID)
    
    if (nrow(mod_nodes) == 0 || nrow(mod_edges) == 0) {
        message(sprintf("Module %s has no eligible nodes or edges; skipping.", mod))
        all_paths[[as.character(mod)]] <- list()
        next
    }
    
    g_mod <- graph_from_data_frame(mod_edges, vertices = mod_nodes, directed = FALSE)
       
    micro_nodes <- V(g_mod)[Group == "Microbiome"]$name
    
    if (length(micro_nodes) == 0) {
        message(sprintf("Module %s has no Microbiome nodes; skipping.", mod))
        all_paths[[as.character(mod)]] <- list()
        next
    }
    
    keep_paths <- list()
    for (start_node in micro_nodes) {
        paths_mod <- all_simple_paths(g_mod, from = start_node, cutoff = 5)
        if (length(paths_mod) == 0) {
            next
        }
        
        paths_named <- lapply(paths_mod, function(p) names(p))
    
        for (p in paths_named) {
            p_group <- mod_nodes$Group[match(p, mod_nodes$ID)]
            
            if (any(p_group == "Microbiome") &&
                any(p_group %in% c("Metabolome", "Lipidome")) &&
                any(p_group == "Proteome")) {
                            
                if (length(p) >= 3) {
                    if (length(p) == 3) {
                        triplets <- matrix(p, nrow = 1)
                    } else {
                        triplets <- embed(rev(p), 3)[, 3:1, drop = FALSE]  
                    }
                    
                    for (row in seq_len(nrow(triplets))) {
                        trip <- triplets[row, ]
                        g_trip <- mod_nodes$Group[match(trip, mod_nodes$ID)]
                        
                        valid <- any(
                            all(g_trip == c("Microbiome", "Metabolome", "Proteome")),
                            all(g_trip == c("Proteome", "Metabolome", "Microbiome")),
                            all(g_trip == c("Microbiome", "Lipidome", "Proteome")),
                            all(g_trip == c("Proteome", "Lipidome", "Microbiome"))
                        )
                        
                        if (valid) {
                            mid <- trip[2]
                            ends <- sort(c(trip[1], trip[3]))
                            path_id <- paste(mid, paste(ends, collapse = "_"), sep = "|")
                            keep_paths[[path_id]] <- trip
                        }
                    }
                }
            }
        }
    }
    if (length(keep_paths) == 0) {
        message(sprintf("Module %s produced no qualifying cross-omics paths.", mod))
    }
    all_paths[[as.character(mod)]] <- keep_paths 
}

# Export data.frame
all_paths <- Filter(function(x) length(x) > 0, all_paths)
if (length(all_paths) == 0) {
  message("No qualifying cross-omics paths were found across all modules; writing an empty Subpath.csv.")
  df_paths1 <- data.frame(Source = character(), Middle = character(), Target = character())
} else {
  final_paths <- unique(unlist(all_paths, recursive = FALSE))
  df_paths1 <- bind_rows(lapply(final_paths, function(x) {
      data.frame(Source = x[1], Middle = x[2], Target = x[3])
  }))
}

write.csv(df_paths1, file.path(outdir, "Subpath.csv"), row.names = FALSE)
