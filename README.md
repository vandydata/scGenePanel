# scGenePanel

[![R-CMD-check](https://github.com/vandydata/scGenePanel/workflows/R-CMD-check/badge.svg)](https://github.com/vandydata/scGenePanel/actions)
[![Bioconductor candidate](https://img.shields.io/badge/Bioconductor-candidate-blue.svg)](https://bioconductor.org/packages/devel/bioc/)


scGenePanel creates a multi-panel visualization for various gene expression metrics on user defined gene, sample group and cell types using single-cell RNAseq data.

## Features



1. A multi-panel plot that contains:
   1. UMAP : Visualize cluster of cells, split by cell identity of interest and highlighted cell type of interest 
   2. Violin plot : Visualize gene expression counts, colored by cell identity of interest and split by cell type of interest 
   3. Tabular plot : Quantify cell counts, cell frequency ratio and quantiles of gene expression counts per cell identity of interest
2. Accepts two commonly used input data object (Seurat, SingleCellExperiment) via `object `parameter
3. "Tableau" or "RColorBrewer" discrete qualitative color palettes options available via `col.palette` parameter
4. Order of cell identity of interest can be changes via `group_order` parameter 
5. Shiny version to explore gene expression data interactively


## Installation

From Github repository (once publicly available)

```R
library(devtools)
install_github("vandydata/scGenePanel")
```

From local Git clone

```R
install.packages("/path/to/scGenePanel", repos = NULL, type = "source")
# with dependencies
install.packages("/path/to/scGenePanel", repos = NULL, type = "source", dependencies = TRUE)
```

## Development

```
# at startup
devtools::load_all() # from package dir, loads but doesn't install

# after edits
devtools::document()
devtools::load_all()
```

