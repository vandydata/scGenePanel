# Preparation record for bundled example data
#
# The files in inst/extdata are intentionally small example objects for the
# vignette, package examples, and unit tests. They are not a complete analysis
# data set and should not be used for biological inference.
#
# Source study:
# Shrestha S, et al. Combinatorial transcription factor profiles predict mature
# and functional human islet alpha and beta cells. JCI Insight. 2021;6(22).
# doi:10.1172/jci.insight.151621. Sequencing data: GEO GSE183568.
#
# Preparation summary:
# 1. A processed healthy human pancreatic-islet Seurat object derived from the
#    source study was used during package development.
# 2. A random 500-cell subset was retained to make examples and tests fast.
# 3. The subset was saved as inst/extdata/human_panc_islets.Rds.
# 4. An equivalent SingleCellExperiment representation was saved as
#    inst/extdata/sce_object.Rds.
#
# The historical random-selection seed was not retained. The distributed files
# are therefore the canonical examples; this script records their provenance
# and preparation rather than reproducing their exact cell identities.
#
# Development intermediates were kept outside the installed package as
# dev/panc.filtered.celltypes.diet.Rds and dev/panc_sub.Rds.
