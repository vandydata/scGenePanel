1. A single "image/plot" that has key information from multiple levels. 
   1. UMAP : Visualize cluster of cells, split by cell identity of interest and highlighted cell type of interest 
   2. Violin plot : Visualize gene expression counts, colored by cell identity of interest and split by cell type of interest 
   3. Tabular plot : Cell counts, cell frequency ratio and quantiles of gene expression counts per cell identity of interest
2. Accepts 2 commonly used input data object (Seurat, SingleCellExperiment) via `object `parameter
3. "Tableu" or "RColorBrewer" discrete qualitative color palettes options available via `col.palette` parameter
4. Order of cell identity of interest can be changes via `group_order` parameter 
5. Shiny version to explore gene expression data interactively