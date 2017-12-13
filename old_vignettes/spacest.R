## ---- eval = FALSE-------------------------------------------------------
#  library(devtools)
#  install_github("ludvigla/ST_analysis3D")

## ---- message=FALSE, warning=FALSE, eval = FALSE-------------------------
#  library(STanalysis3D)
#  data(bcST)

## ---- echo=FALSE, results='asis'-----------------------------------------
# Check the format of te first expression dataset
#knitr::kable(head(bcST[[1]][, 1:8]))

## ---- eval = FALSE-------------------------------------------------------
#  # Replace ensembl gene ids with HGNC symbols
#  exp.list <- ensembl2hgnc(exp.list)
#

## ---- eval = FALSE-------------------------------------------------------
#  # Initialize spaceST object with batch correction
#  spaceST.object <- get.spaceST(bcST, correction = T)
#  spaceST.object

## ---- fig.show='hold', eval = FALSE--------------------------------------
#  # QC plot
#  plot.QC.spaceST(spaceST.object)
#
#  # PCA plot
#  pca.spaceST(spaceST.object)

## ---- fig.show='hold', fig.width=6, eval = FALSE-------------------------
#  # Or if you want to split up relicates in separate plots
#  plot.QC.spaceST(spaceST.object, separate = T)

## ---- eval = FALSE-------------------------------------------------------
#  ST.topics <- topic.compute(spaceST.object)

## ---- eval=FALSE, echo=FALSE, results='asis'-----------------------------
#  # Check the format of te first expression dataset
#  knitr::kable(head(ST.topics[, 1:8]))

## ---- echo=FALSE, results='asis'-----------------------------------------
knitr::kable(head(mtcars, 10))

