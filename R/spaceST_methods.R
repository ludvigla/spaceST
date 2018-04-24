#' @import methods
NULL

#' The spaceST class.
#'
#' @description An S4 class object used to filter and store spatial transcriptomics data from
#' multiple tissue sections
#' @slot expr Object of class dcGMatrix containing gene expression values.
#' @slot norm.data Object of class dcGMatrix containing normalized expression values.
#' @slot lda.results Object of class topics obtained using the \link[cellTree]{compute.lda} function from cellTree.
#' @slot tsne Matrix object representation of t-SNE results.
#' @slot snn Matrix object representation of the SNN graph.
#' @slot coordinates Data.frame representation of section ID and array coordinates.
#' @slot filter.settings List of filtering settings.
#' @slot meta.data List of meta data.
#' @slot version package version.
#' @name spaceST
#' @rdname spaceST
#' @aliases spaceST-class
#' @exportClass spaceST
spaceST <- setClass(
  "spaceST",
  slots = c(
    expr = "dgCMatrix",
    norm.data = "dgCMatrix",
    lda.results = "ANY",
    reducedDims = "ANY",
    tsne = "matrix",
    snn = "ANY",
    coordinates = "data.frame",
    filter.settings = "list",
    meta.data = "list",
    version = "ANY",
    status.expr = "character"
  )
)


#' Show method for spaceST.
#'
#' @param object Object of class spaceST.
#' @aliases spaceST-show
#' @name show
#' @rdname show-methods
setMethod("show", signature = "spaceST", definition = function(object) {
  cat("An object of class ", class(object), "\n\n", sep = "")
  cat("Filter settings:\n")
  cat("  Filtered features with less than ", object@filter.settings[[1]], " unique genes \n", sep = "")
  cat("  Filtered genes with expression values =< ", object@filter.settings[[2]], " in =< ", object@filter.settings[[3]],
      " features\n\n", sep = "")
  cat("expr status:", object@status.expr)
})


#' Batch correction of spaceST data.
#'
#' Function used to run a batch correction between samples.
#' @description Run batch correction of replicate sections.
#' @param object Object of class spaceST.
#' @seealso \link[countClust]{BatchCorrectedCounts}
#' @importFrom CountClust BatchCorrectedCounts
#' @return Batch corrected data frame in slot corrected
#' @rdname batch.correct
#' @export
setGeneric(
  name = "batch.correct",
  def = function(object) standardGeneric("batch.correct")
)
#' @rdname batch.correct
#' @export
setMethod(
  f = "batch.correct",
  signature = "spaceST",
  definition = function(object) {
    if (class(object) != "spaceST"){
      stop("Wrong input format.")
    }
    object@expr <- as(t(BatchCorrectedCounts(t(as.matrix(object@expr)), object@coordinates[, 1], use_parallel = T)), "dgCMatrix")
    object@status.expr <- "batch corrected"
    return(object)
  }
)


#' Normalize ST data.
#'
#' Normalization of gene expression data using either counts per 10k or scran.
#' @param object Object of class spaceST.
#' @param method Select normalization method, default "cp10k" [options: "cp10k", "scran"].
#' @param log2 Logical specifying whether or not data should be log2-transformed.
#' @param clusters Integer vector specifying clusters used for normalization with scran on large gene expression matrices.
#' By default, a clustering step is performed using \link[scran]{quickCluster} if there are more than 500 array spots in the
#' expression matrix.
#' @seealso \link[scater]{Normalize} and \link[scran]{computeSumFactors}
#' @importFrom scran quickCluster computeSumFactors
#' @rdname normalize
#' @return Normalized expression data.
#' @export
setGeneric(name = "NormalizespaceST",
           def = function(object,
                          method = "scran",
                          log2 = F,
                          clusters = NULL, ...
           ) standardGeneric("NormalizespaceST"))
#' @rdname normalize
#' @export
setMethod(
  f = "NormalizespaceST",
  signature = "spaceST",
  definition = function(object,
                        method = "cp10k",
                        log2 = F,
                        pcount = 1,
                        cluster = NULL, ...) {
    if (class(object) != "spaceST"){
      stop("Wrong input format.")
    }
    expr.data <- as.matrix(object@expr)
    if (log2 & method == "cp10k") {
      object@norm.data <- log2(calc_cpm(expr.data) + 1)
    } else {
      if (method == "scran") {
        sce <- SingleCellExperiment(assays = list(counts = expr.data, logcounts = log2(expr.data + pcount)))
        if (ncol(object@expr) > 500 | !is.null(cluster)) {
          if (!is.null(cluster)) {
            sce <- computeSumFactors(sce, cluster = cluster)
          } else {
            q.clust <- quickCluster(as.matrix(expr.data))
            sce <- computeSumFactors(sce, cluster = q.clust)
          }

          object@meta.data$size_factor <- sizeFactors(sce)
        } else {
          sce <- computeSumFactors(sce)
          object@meta.data$size_factor <- sizeFactors(sce)
        }
        if (sum(sizeFactors(sce) == 0) > 0) {
          stop("Zero sum factors not allowed. Filter data before normalization with scran.")
        }
        sce <- normalize(sce)
        object@norm.data <- as(logcounts(sce), "dgCMatrix")
      } else if (method == "cp10k") {
        object@norm.data <- as(calc_cp10k(as.matrix(expr.data)), "dgCMatrix")
      }
    }
    return(object)
  }
)


#' Plot unique genes per feature and transcripts per feature.
#'
#' his function is used to plot the unique genes per feature and transcipts per feature
#' distributions as histograms. By default, the function will take the batch corrected dataset
#' if present.
#' @param object Object of class spaceST.
#' @param separate Logical indicating whether replicates should be ploteed separ2ately
#' @rdname QC
#' @return Histograms of unique genes per feature and transcripts per feature distributions.
#' @export
setGeneric("plot.QC.spaceST", function(object, separate = F) standardGeneric("plot.QC.spaceST"))
#' @rdname QC
#' @docType methods
#' @export
setMethod("plot.QC.spaceST", "spaceST", function(object, separate = F) {
  if (class(object) != "spaceST") {
    stop("Wrong input format")
  }
  if (length(object@expr) > 0) {
    df <- ST_statistics(as.matrix(object@expr))
    ST_statistics_plot(df, separate)
  }
})


#' PCA method.
#'
#' Run PC analysis on spaceST object.
#' @param object Object of class spaceST.
#' @param ntop Number of variable genes to select for PC analysis.
#' @param ncomponents Number of components that should be returned.
#' @param exprs_values Select target object to use as input for PC analysis [options: "norm.data", "expr"]
#' @param ... Parameters passed to \link[stats]{prcomp}.
#' @seealso \link[stats]{prcomp}
#' @aliases spaceST-pca
#' @export
setGeneric("spPCA", function(object, ntop = 500, ncomponents = 2, exprs_values = "norm.data", ...) standardGeneric("spPCA"))
#' @export
setMethod(
  "spPCA",
  signature = "spaceST",
  definition = function(
  object,
  ntop = 500,
  exprs_values = "norm.data",
  ...
) {
  if (exprs_values == "norm.data") {
    input.data <- as.matrix(object@norm.data)
  } else if (exprs_values == "exprs") {
    input.data <- as.matrix(object@expr)
  }
  high.genes <- names(sort(apply(input.data, 1, var), decreasing = T)[1:ntop])
  input.data <- input.data[high.genes, ]
  pca <- prcomp(t(input.data), scale. = T, center = T, ...)
  pcs <- pca$x
  object@reducedDims <- pcs
  return(object)
})
