library(scGenePanel)
library(Seurat)


# make changes
# Click on Build > Install
# This will clean up/remove existing install and re-install, takes a few seconds
# and then reload session files.
#
# Use this to load functions from the package so
# they can be tested outside of rebuilding package
# devtools::load_all()


# Load the data
#panc <- readRDS("panc.filtered.celltypes.diet.Rds")
#panc_sub <- panc[, sample(colnames(panc), size =3000, replace=F)]
#rm(panc)
#saveRDS(panc_sub, "panc_sub.Rds", compress=FALSE)
panc_sub <- readRDS("dev/panc_sub.Rds")

library(scGenePanel)
create_gene_panel(panc_sub,
                   "INS",
                   "status",
                   "Beta",
                   "CellTypes")






devtools::load_all()
umap_panel(panc_sub,
           cell_type_colname="CellTypes",
           cell_type_name="Beta",
           meta_group="sex",
           gene="GCG",
           levels_idents,
           group_order=c("Female", "Male")
           )

umap_panel(panc_sub,
           cell_type_colname="CellTypes",
           cell_type_name="Beta",
           meta_group="sex",
           gene="GCG",
           levels_idents)


violin_panel(panc_sub,
             cell_type_colname="CellTypes",
             cell_type_name="Beta",
             meta_group="sex",
             gene="GCG",
             col_palette="Tableau",
             levels_idents)

boxplot_panel(panc_sub,
             cell_type_colname="CellTypes",
             cell_type_name="Beta",
             meta_group="sex",
             gene="GCG",
             col_palette="Tableau",
             levels_idents)

cellfreq_panel(panc_sub,
               cell_type_colname="CellTypes",
               cell_type_name="Beta",
               meta_group="sex",
               gene="GCG",
               col_palette="Tableau",
               levels_idents)


cellfreq2_panel(panc_sub,
               cell_type_colname="CellTypes",
               cell_type_name="Beta",
               meta_group="sex",
               gene="GCG",
               col_palette="Tableau",
               levels_idents)

## Small test

panc_human_test <- readRDS("D:/data/CDS-datasets/scRNA-Seq/human_pancreas-islets_seurat.Rds")  # Great for testing, small and good cluster separation

object=panc_human_test
cell_type_colname="CellTypes"
cell_type_name="Alpha"
meta_group="Age"
group_order=c("14y", "39y", "50y", "59y", "66y")
gene="GCG"
col_palette = "Tableau"
levels_idents <- c("39y", "50y", "14y", "59y", "66y")
seurat_obj=object
Seurat::Idents(seurat_obj) <- meta_group


library(scGenePanel)
create_gene_panel(panc_human_test,
                  "GCG",
                  "Age",
                  "Alpha",
                  "CellTypes")

devtools::load_all()
umap <- umap_panel(object,
           cell_type_colname="CellTypes",
           cell_type_name="Alpha",
           meta_group="Age",
           gene="GCG",
           levels_idents)

violin <- violin_sig_panel(object,
                 cell_type_colname="CellTypes",
                 cell_type_name="Alpha",
                 meta_group="Age",
                 gene="GCG",
                 col_palette="Tableau",
                 levels_idents)

table <- cellfreq_panel(object,
                cell_type_colname="CellTypes",
                cell_type_name="Alpha",
                meta_group="Age",
                gene="GCG",
                col_palette="Tableau"
)

cellfreq2_panel(object,
                cell_type_colname="CellTypes",
                cell_type_name="Alpha",
                meta_group="Age",
                gene="GCG",
                col_palette="Tableau"
)

object=panc_human_test
cell_type_colname="CellTypes"
cell_type_name="Alpha"
meta_group="Age"
group_order=c("Female", "Male")
gene="INS"
col_palette = "Tableau"
levels_idents <- c("Male", "Female")
levels_idents <- c("39y", "50y", "14y", "59y", "66y")
seurat_obj=object
Seurat::Idents(seurat_obj) <- meta_group


# to rebuild docs
devtools::document()



