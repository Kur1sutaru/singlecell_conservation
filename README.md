<p align="center">
  <img src="man/figures/logo.png" alt="singleCellConservation logo" width="180"/>
</p>

<h1 align="center">singleCellConservation</h1>

<p align="center">
  <img src="https://img.shields.io/badge/R%20package-in%20development-blue" alt="status"/>
  <img src="https://img.shields.io/badge/license-MIT-green" alt="license"/>
</p>

Cross-species conservation analysis for single-cell marker genes, built on
Seurat. Given a Seurat object with a cell-subtype column and a species
column, this package runs a 4-step pipeline:

1. **Unfiltered `FindConservedMarkers()`** across all subtypes (meta-analysis
   across species, no thresholds applied yet).
2. **Classification** of each gene into `Conserved` (passes independently in
   *every* species), `Partially conserved` (passes in 2+ but not all),
   `Species-specific` (passes in exactly 1), or `Not significant` (passes in
   none), using per-species log2FC, p-value (raw or adjusted), and detection
   rate (`pct.1`) thresholds.
3. **Pairwise Pearson correlation** of mean expression between every pair of
   species, per subtype.
4. **Expression stability** scoring: normalized SD of expression rank across
   species, bucketed into Highly stable / Moderately stable / Variable /
   Highly variable.

This was originally developed for a cross-species tongue atlas
(Human/Marmoset/Mouse) but works on any Seurat object with the right
metadata columns.

## Installation

```r
# install.packages("devtools")
devtools::install_github("Kur1sutaru/singlecell_conservation")
```

## Usage

```r
library(singleCellConservation)
library(Seurat)

results <- run_conservation_pipeline(
  obj             = my_seurat_object,
  subtype_col     = "cca_subtypes",
  species_col     = "species",
  species_list    = c("Human", "Marmoset", "Mouse"),
  logfc_thresh    = 0.6,
  pval_thresh     = 0.05,
  min_pct         = 0.25,
  pval_col_suffix = "_p_val_adj"   # or "_p_val" for raw p-values
)

# results$step1_raw          -- unfiltered FindConservedMarkers output
# results$step2_classified   -- every gene classified, with a "conservation" column
# results$step3_pearson      -- pairwise Pearson r per subtype
# results$step4_stability    -- stability classification of conserved genes

conserved_strict <- results$step2_classified |>
  dplyr::filter(conservation == "Conserved")

plot_conserved_count(conserved_strict)
plot_pearson_heatmap(results$step3_pearson)
plot_stability(results$step4_stability)
```

### Renaming subtype labels before running the pipeline

```r
my_seurat_object <- rename_subtypes(
  my_seurat_object,
  subtype_col = "cca_subtypes",
  rename_map = c(
    "NK/Cytotoxic T cells" = "CD8+ Cytotoxic T cells",
    "\u03b3\u03b4 T cells" = "Th17-like T cells",
    "AXL+ Dendritic cells" = "Dendritic cells"
  )
)
```

## Individual pipeline steps

Each step can also be run independently:

```r
step1 <- run_step1_conserved_markers(obj, "cca_subtypes", "species", c("Human","Marmoset","Mouse"))
step2 <- classify_step2(step1, c("Human","Marmoset","Mouse"), logfc_thresh = 0.6,
                         pval_thresh = 0.05, pval_col_suffix = "_p_val_adj", min_pct = 0.25)
step3 <- compute_pearson_pairwise(obj, step2[step2$conservation == "Conserved", ],
                                   c("Human","Marmoset","Mouse"), "cca_subtypes", "species")
step4 <- compute_stability(step2[step2$conservation == "Conserved", ], attr(step2, "fc_cols"))
```

## A note on `total_step1` reconciliation

Because Step 1 is unfiltered, `total_step1` (all genes tested for a subtype)
will generally be **larger** than `Conserved + Partially conserved +
Species-specific`, since some genes pass in zero species and land in
`Not significant`. If you're building your own summary tables downstream,
make sure to include that fourth bucket or your totals won't reconcile.

## Development

```r
devtools::load_all()
devtools::test()
devtools::document()   # regenerates NAMESPACE / man/ from roxygen comments
```

## License

MIT
