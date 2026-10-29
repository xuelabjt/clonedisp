# Generate CloneDisp labels from inferCNV

CloneDisp does **not** call CNVs and does not choose clones. It consumes a
cell-level table:

```text
sample,state,clone
T1,e1,C1
T1,e1,C2
```

The normal workflow is:

```text
raw RNA counts + normal references
        |
        v
R inferCNV denoised CNV matrix
        |
        v
PCA and clustering within each patient
        |
        v
sample,state,clone labels
        |
        v
CloneDisp
```

The transcriptomic `state` can come directly from the single-cell object. The
`clone` label must come from clustering the CNV profile, not from clustering
RNA expression.

## 1. Run R inferCNV

Pass **raw counts** to inferCNV. Depth normalization, log transformation,
reference centering, smoothing, and denoising happen inside
`infercnv::run()`. Do not pass a log-normalized RNA matrix.

Minimal R structure:

```r
library(infercnv)

cnv <- CreateInfercnvObject(
  raw_counts_matrix = counts,
  gene_order_file = gene_order,
  annotations_file = annotations,
  ref_group_names = c("reference"),
  chr_exclude = c("chrX", "chrY", "chrM")
)

cnv <- infercnv::run(
  cnv,
  cutoff = 0.05,
  out_dir = "infercnv_run",
  window_length = 101,
  smooth_method = "pyramidinal",
  cluster_by_groups = FALSE,
  max_centered_threshold = 3,
  scale_data = FALSE,
  HMM = FALSE,
  analysis_mode = "samples",
  denoise = TRUE,
  save_final_rds = TRUE
)
```

`cutoff` is data-dependent. The value `0.05` is the setting used for this
sparse 10x PDAC project; the inferCNV default of `1` removed too many genes.
Always check that the final object retains the genes required for your
analysis.

## 2. Cluster the denoised CNV matrix per patient

Prepare a metadata table with one row per tumor cell:

```text
cell,sample,state
T1_AAC...,T1,e1
```

Then run the bundled R helper:

```bash
Rscript scripts/infercnv_to_labels.R \
  infercnv_run/run.final.infercnv_obj \
  cell_metadata.csv \
  clonedisp_labels \
  0.5,0.6,0.7,0.8,0.9,1.0 \
  5 30 20
```

The arguments are:

```text
FINAL_INFERCNV_RDS METADATA_CSV OUTPUT_DIR
RESOLUTIONS MIN_CLONE_CELLS N_PCS N_NEIGHBORS
```

The helper:

1. takes the inferCNV `expr.data`;
2. converts the reference-centered fold changes with `log2`;
3. removes zero-variance genes;
4. performs PCA and Leiden clustering separately in each patient;
5. reassigns microclusters below `MIN_CLONE_CELLS` to the most correlated
   retained clone;
6. writes one `sample,state,clone,cell` table per resolution;
7. writes a `definitions.csv` manifest for `clonedisp-multi`.

The exact clustering choices are analysis decisions and should be reported.
CloneDisp tests the metric across definitions; it does not make an unsuitable
clone definition valid.

## 3. Run CloneDisp across resolutions

```bash
clonedisp-multi \
  --manifest clonedisp_labels/definitions.csv \
  --output-prefix clonedisp_results \
  --min-clone-cells 5 \
  --min-state-cells 20 \
  --plot
```

## infercnvpy

The same input rules apply to Python `infercnvpy`: run
`infercnvpy.tl.infercnv` to obtain a CNV profile, cluster those CNV profiles in
each patient, and export `sample,state,clone`. R inferCNV and infercnvpy labels
must remain separate analyses; their raw clone labels should not be mixed.

## Common mistakes

- Clustering RNA expression and calling the clusters CNV clones.
- Using normalized expression instead of raw counts as inferCNV input.
- Mixing normal/reference cells into the clone-clustering set.
- Pooling patients before clone calling when clones are intended to be
  patient-specific.
- Treating a clone number or resolution as a natural biological constant.
- Combining labels from multiple resolutions into one `K_eff` calculation.
