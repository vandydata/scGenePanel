# scGenePanel example data

`human_panc_islets.Rds` and `sce_object.Rds` are equivalent, small example
representations of a processed healthy human pancreatic-islet single-cell
RNA-seq data set. They are included only for package examples, the vignette,
and automated tests.

## Provenance

The examples were derived from the healthy human-islet scRNA-seq data reported
by Shrestha S, et al. *Combinatorial transcription factor profiles predict
mature and functional human islet alpha and beta cells*. JCI Insight.
2021;6(22). doi:10.1172/jci.insight.151621. Sequencing data are available
through GEO accession GSE183568.

During package development, a random 500-cell subset of a processed object
derived from this study was retained to keep package examples and tests fast.
The historic random-selection seed was not retained; the distributed files are
the canonical example objects. They must not be used for biological inference.
See `inst/scripts/prepare-example-data.R` for the preparation record.

## Files

| File | Class | Dimensions | Assays and reductions |
|---|---|---:|---|
| `human_panc_islets.Rds` | `Seurat` | 26,107 genes x 500 cells | `RNA`; PCA and UMAP |
| `sce_object.Rds` | `SingleCellExperiment` | 26,107 genes x 500 cells | `counts`; PCA and UMAP |

Both objects include the following cell-level metadata columns:

`orig.ident`, `nCount_RNA`, `nFeature_RNA`, `Library_ID`, `Age`, `Source`,
`Sex`, `Batch`, `Replicates`, `percent.mt`, `RNA_snn_res.0.6`,
`seurat_clusters`, and `CellTypes`. The Seurat object also includes
`RNA_snn_res.0.5`.
