## ---- eval = FALSE-------------------------------------------------------
#  library(devtools)
#  install_github("ludvigla/spaceST")

## ---- message=FALSE, warning=FALSE---------------------------------------
library(spaceST)
load(file = system.file("data/MOB.RData", package = "spaceST"))

## ---- echo=FALSE, results='asis'-----------------------------------------
knitr::kable(MOB1[1:10, 1:5])

## ----create spaceST object-----------------------------------------------
# Initialize spaceST object and filter
spST <- CreatespaceSTobject(list(MOB1, MOB2), unique.genes = 300, min.exp = 2, min.features = 10, delimiter = "x", filter.genes = "^Mt-*|Malat1")
spST

dim(spST)

## ---- fig.show='hold', fig.width = 8, fig.height = 5---------------------
# QC plot
plot_QC_spaceST(spST)

## ---- fig.show='hold', fig.width = 11, fig.height = 7--------------------
# Or if you want to split up relicates in separate plots
plot_QC_spaceST(spST, separate = T)

## ---- eval = FALSE-------------------------------------------------------
#  # Initialize spaceST object and filter
#  selection.files <- list("path1", "path2")
#  
#  # Don't forget to specify the delimiter
#  spST.subset <- spots.under.tissue(spST, selection.files, delimiter = "x")
#  spST
#  
#  dim(spST)

## ----gene viz, fig.width = 11, fig.height = 5----------------------------
spatial.heatmap(spST, value = "Nrgn", invert.heatmap = T, size = 3)

## ----gene viz HE, fig.width = 11, fig.height = 5-------------------------
HE.list <- list(system.file("data/HE_rep1.jpg", package = "spaceST"),
                system.file("data/HE_rep2.jpg", package = "spaceST"))
spatial.heatmap(spST, value = "Nrgn", invert.heatmap = T, HE.list = HE.list, alpha = 0.6, size = 3)

## ----normalize, eval = FALSE---------------------------------------------
#  # Normalize using scran
#  spST <- NormalizespaceST(spST, method = "scran")

## ----load, echo=FALSE, results='hide',message=FALSE----------------------
load("~/breast_cancer/spaceST/spST")

## ----gene viz norm, fig.width = 11, fig.height = 5-----------------------
spatial.heatmap(spST, value = "Nrgn", type = "norm.data", invert.heatmap = T, HE.list = HE.list, alpha = 0.6, size = 3)

## ----pca, fig.width = 4, fig.height = 3----------------------------------
spST <- spPCA(spST)
plotPCA(spST, components = c(1, 2))

## ----fig.width = 11, fig.height = 5--------------------------------------
spatial.heatmap(spST, value = "PC1", type = "pca", invert.heatmap = T, HE.list = HE.list, alpha = 0.6, size = 3)

## ---- eval = FALSE-------------------------------------------------------
#  spST <- topic_compute(spST)

## ----fig.width = 11, fig.height = 5--------------------------------------
spatial.heatmap(spST, value = "1", type = "topics", invert.heatmap = T, HE.list = HE.list, size = 3, alpha = 0.6)
spatial.heatmap(spST, value = "2", type = "topics", invert.heatmap = T, HE.list = HE.list, size = 3, alpha = 0.6)
spatial.heatmap(spST, value = "3", type = "topics", invert.heatmap = T, HE.list = HE.list, size = 3, alpha = 0.6)

## ------------------------------------------------------------------------
top_features <- ExtractTopFeaturesST(spST)
knitr::kable(top_features[1:10, ])

## ----clusters on tissue, fig.width = 11, fig.height = 5------------------
spatial.clusters(spST, HE.list = HE.list, size = 3, alpha = 0.6)

