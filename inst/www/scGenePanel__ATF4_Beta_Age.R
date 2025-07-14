# Data source: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE183568
data <- subset(x = data, idents = c("Beta", "Alpha",  "Ductal", "Acinar")) #  10,000 cells

# Create figure panel
create_gene_panel(data,
                  gene="ATF4",
                  meta_group = "Age",
                  cell_type_name = "Beta",
                  cell_type_colname = "CellTypes")
