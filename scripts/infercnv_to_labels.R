#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(irlba))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop(
    paste(
      "Usage:",
      "Rscript infercnv_to_labels.R FINAL_INFERCNV_RDS METADATA_CSV OUTPUT_DIR",
      "[RESOLUTIONS] [MIN_CLONE_CELLS] [N_PCS] [N_NEIGHBORS]"
    )
  )
}

infercnv_path <- args[[1]]
metadata_path <- args[[2]]
output_dir <- args[[3]]
resolutions <- if (length(args) >= 4) {
  as.numeric(strsplit(args[[4]], ",", fixed = TRUE)[[1]])
} else {
  c(0.5, 0.7, 1.0)
}
min_clone_cells <- if (length(args) >= 5) as.integer(args[[5]]) else 5L
n_pcs <- if (length(args) >= 6) as.integer(args[[6]]) else 30L
n_neighbors <- if (length(args) >= 7) as.integer(args[[7]]) else 20L

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cnv_object <- readRDS(infercnv_path)
if (!("expr.data" %in% slotNames(cnv_object))) {
  stop("The RDS object does not expose an inferCNV expr.data slot")
}

metadata <- read.csv(
  metadata_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
required_columns <- c("cell", "sample", "state")
missing_columns <- setdiff(required_columns, colnames(metadata))
if (length(missing_columns)) {
  stop(
    "Metadata is missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}
if (anyDuplicated(metadata$cell)) {
  stop("Metadata cell identifiers must be unique")
}
missing_cells <- setdiff(metadata$cell, colnames(cnv_object@expr.data))
if (length(missing_cells)) {
  stop(length(missing_cells), " metadata cells are absent from inferCNV expr.data")
}

refine_small_clusters <- function(expr, labels, min_cells) {
  labels <- as.character(labels)
  repeat {
    sizes <- table(labels)
    small <- names(sizes)[sizes < min_cells]
    large <- names(sizes)[sizes >= min_cells]
    if (!length(small) || length(large) < 2) {
      break
    }
    centers <- sapply(large, function(cluster) {
      Matrix::rowMeans(expr[, labels == cluster, drop = FALSE])
    })
    centers <- as.matrix(centers)
    for (cluster in small) {
      for (cell in which(labels == cluster)) {
        correlations <- apply(
          centers,
          2,
          function(center) cor(expr[, cell], center)
        )
        labels[[cell]] <- names(which.max(correlations))
      }
    }
  }
  labels
}

all_rows <- list()
manifest_rows <- list()

for (sample_name in sort(unique(metadata$sample))) {
  sample_meta <- metadata[metadata$sample == sample_name, , drop = FALSE]
  if (nrow(sample_meta) < 6) {
    warning("Skipping ", sample_name, ": fewer than six tumor cells")
    next
  }

  expression <- log2(
    cnv_object@expr.data[, sample_meta$cell, drop = FALSE]
  )
  expression <- expression[apply(expression, 1, sd) > 0, , drop = FALSE]
  pc_count <- min(n_pcs, nrow(sample_meta) - 1, nrow(expression) - 1)
  pca <- irlba::prcomp_irlba(
    t(expression),
    n = pc_count,
    center = TRUE,
    scale. = FALSE
  )
  rownames(pca$x) <- sample_meta$cell

  dummy <- matrix(
    1,
    nrow = 1,
    ncol = nrow(sample_meta),
    dimnames = list("dummy", sample_meta$cell)
  )
  options(Seurat.object.assay.version = "v3")
  seurat_object <- CreateSeuratObject(counts = Matrix(dummy, sparse = TRUE))
  seurat_object[["cnv_pca"]] <- CreateDimReducObject(
    embeddings = pca$x,
    key = "cnv_",
    assay = DefaultAssay(seurat_object)
  )
  seurat_object <- FindNeighbors(
    seurat_object,
    reduction = "cnv_pca",
    dims = seq_len(ncol(pca$x)),
    k.param = min(n_neighbors, nrow(sample_meta) - 1),
    verbose = FALSE
  )

  for (resolution in resolutions) {
    set.seed(1)
    resolution_key <- gsub("\\.", "p", format(resolution, trim = TRUE))
    cluster_column <- paste0("res_", resolution_key)
    seurat_object <- FindClusters(
      seurat_object,
      resolution = resolution,
      algorithm = 4,
      random.seed = 1,
      cluster.name = cluster_column,
      group.singletons = FALSE,
      verbose = FALSE
    )
    labels <- refine_small_clusters(
      expression,
      seurat_object@meta.data[[cluster_column]],
      min_clone_cells
    )
    labels <- paste0("C", as.integer(factor(labels)))

    table <- data.frame(
      sample = sample_name,
      state = sample_meta$state,
      clone = labels,
      cell = sample_meta$cell,
      stringsAsFactors = FALSE
    )
    definition <- paste0("res ", format(resolution, trim = TRUE))
    filename <- paste0("labels_", cluster_column, ".csv.gz")
    write.csv(
      table,
      gzfile(file.path(output_dir, filename)),
      row.names = FALSE
    )
    all_rows[[length(all_rows) + 1L]] <- cbind(
      resolution = format(resolution, trim = TRUE),
      table
    )
    manifest_rows[[length(manifest_rows) + 1L]] <- data.frame(
      definition = definition,
      path = filename,
      stringsAsFactors = FALSE
    )
  }
}

if (!length(all_rows)) {
  stop("No sample produced labels")
}

write.csv(
  do.call(rbind, all_rows),
  gzfile(file.path(output_dir, "labels_all_resolutions.csv.gz")),
  row.names = FALSE
)
write.csv(
  do.call(rbind, manifest_rows),
  file.path(output_dir, "definitions.csv"),
  row.names = FALSE
)

cat("Wrote CloneDisp labels and definitions.csv to", output_dir, "\n")
