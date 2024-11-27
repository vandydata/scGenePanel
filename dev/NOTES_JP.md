## Development - in `.gitignore`





### Must-do checks

1. R CMD build
   ```
   cd d:\data\scGenePanel
   R CMD build .
   ```

2. Check with codetools (I think same results as lintr)
   ```
   library(scGenePanel)
   library(scGenePanel)
   BiocCheck::BiocCheckGitClone()
   
   checkUsagePackage("scGenePanel", all=TRUE)
   ```



### Bioconductor Checks

https://contributions.bioconductor.org/general.html

Tested following on 2023-07-25 and all passed.

```
# 1. Build it without errors or warnings
cd cd d:\data\scGenePanel
R CMD build .

# 2. BiocCheck
library(scGenePanel)
checkUsagePackage("scGenePanel", all=TRUE)
BiocCheck::BiocCheckGitClone()

# 3. If any ERROR, ARNINGS and NOTES left, we need to justify them

# 4. File names - Do not use filenames that differ only in case, as not all file systems are case-sensitive.

# 5. source package from BUILD should be <5MB on disk

# 6. Build time <10 mins

# 7. Memory < 3GB (tested in 32-bit R 4.1.2)

# 8. Individual file sizes < 5MB

# 9. Exclude undesirable files via .gitignore and .Rbuildignore
```

Knit vignette from Rmd to HTML, and commit it.













.Rproj.user
*.Rproj
*.Rhistory

While private, one can provide GitHub PAT and install it.

```
library(devtools)
install_github("vandydata/scGenePanel", auth_token="****************************")

```

From local Git clone

```R
install.packages("D:/data/scGenePanel", repos = NULL, type = "source")
# with dependencies
install.packages("D:/data/scGenePanel", repos = NULL, type = "source", dependencies = TRUE)
```

### Reload files without re-installing package

This will reload all source files in the package directly from disk

```
devtools::load_all("/path/to/scGenePanel/")
devtools::check("/path/to/scGenePanel/")

devtools::load_all("d:/data/scGenePanel")
devtools::check("d:/data/scGenePanel")
```

Reloads and reruns roxygen documentation, runs R CMD build, and reloads the package.

```
devtools::document("d:/data/scGenePanel/")
```

final test:

```
R CMD build .
remove.packages("scGenePanel")  
R CMD INSTALL mypackage_1.0.tar.gz
```





### Check:

```
library(BiocCheck)
#BiocCheck("<packageDirOrTarball>")
BiocCheck("")
```

In Rstudio

```sh
lintr:lint("path/to"file)
```

## check with codetools

```
library(scGenePanel)
checkUsagePackage("scGenePanel", all=TRUE)
```

output looks like:

```
cellfreq_panel : <anonymous>: parameter 'e' may not be used
cellfreq_panel : <anonymous>: local variable 'selected_cells_expressed' assigned but may not be used
cellfreq_panel: no visible global function definition for 'unit'
cellfreq_panel: parameter 'col_palette' may not be used
create_gene_panel : <anonymous>: parameter 'e' may not be used
create_gene_panel : <anonymous>: local variable 'selected_cells_expressed' assigned but may not be used
discrete_col_palette: local variable 'pal' used as function with no apparent local function definition
make_seurat : <anonymous>: parameter 'e' may not be used
umap_panel: parameter 'col_palette' may not be used
violin_panel : <anonymous>: parameter 'e' may not be used
violin_panel : <anonymous>: local variable 'selected_cells_expressed' assigned but may not be used
violin_panel: local variable 'panel_figure' assigned but may not be used
```



### Build docs

```
roxygen2::roxygenise() 
```

