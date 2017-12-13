#' @title The spaceST class
#'
#' @description An S4 class object used to filter and store spatial transcriptomics data from
#' multiple tissue sections
#' @slot raw.data Raw input data
#' @slot expr A data.frame with filtered gene expression data
#' @slot corrected A data.frame with batch corrected gene expressiond data
#' @slot normalized A matrix with normalized expression values
#' @slot topics Topic matrix obtained using compute.lda() function from cellTree
#' @slot tsne Tsne results in matrix format
#' @slot snn SNN adjacency matrix
#' @slot coordinates A data.frame to represent feature coordinates in the expression dataset
#' @slot filter.settings A list with settings used for filtering
#' @slot meta.data Meta data slot
#' @slot version package version
#' @name spaceST
#' @rdname spaceST
#' @aliases spaceST-class
#' @exportClass spaceST
spaceST <- setClass(
  "spaceST",
  slots = c(
    raw.data = "ANY",
    expr = "data.frame",
    corrected = "data.frame",
    normalized = "matrix",
    topics = "matrix",
    tsne = "matrix",
    snn = "ANY",
    coordinates = "data.frame",
    filter.settings = "list",
    meta.data = "list",
    version = "ANY"
  )
)

#' show method for spaceST
#'
#' @param object A spaceST object
#' @aliases spaceST-show
#' @name show
#' @rdname show-methods
setMethod("show", signature = "spaceST", definition = function(object) {
  cat("An object of class ", class(object), "\n\n", sep = "")
  cat("Filter settings:\n")
  cat("  Filtered features with less than ", object@filter.settings[[1]], " unique genes \n", sep = "")
  cat("  Filtered genes with expression values =< ", object@filter.settings[[2]], " in =< ", object@filter.settings[[3]],
      " features\n\n", sep = "")
  cat("Summary statistics of")
  if (length(object@corrected) > 0) {
    cat(" corrected data: \n")
    print(summary(ST_statistics(object@corrected[, 1:2])[, 1:2]))
  } else {
    cat(" filtered data: \n")
    print(summary(ST_statistics(object@expr[, 1:2])[, 1:2]))
  }
  if (length(object@topics) > 0) {
    cat("\nLDA results:\n")
    cat("  Selected ", ncol(object@topics), " topics from a ", ncol(object@expr), " 'document' collection\n", sep = "")
    tab <- table(object@meta.data$clusters)
    tab <- paste("  ", names(tab), tab, sep = "\t\t")
    cat("\tcluster:\tnumber of features:\n")
    cat(tab, sep = "\n")
    cat("  Distance method used for clustering:\t", object@meta.data$dist.method, "\n", sep = "")
    cat("  Tree method used for clustering:\t", object@meta.data$tree.method, "\n", sep = "")
    cat("  Minimum cluster size:\t\t\t", object@meta.data$minClusterSize, sep = "")
  }
})

#' Batch correction of spaceST data
#'
#' @param object spaceST object with expression data to correct
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
    object@corrected <- as.data.frame(t(BatchCorrectedCounts(t(object@expr), object@coordinates[, 1], use_parallel = T)))
    return(object)
  }
)

#' Normalize ST data
#'
#' @param object spaceST object with expression data to correct
#' @param method select normalization method, default "cp10k" (options: "cp10k", "scran")
#' @param log2 Logical specifying whether or not data should be log2-transformed
#' @param clusters Integer vector specifying clusters used for normalization of large
#' sce objects. Only for "scran" normalization.
#' @param ... Argumants passed to
#' @importFrom scater newSCESet sizeFactors
#' @importFrom scran quickCluster computeSumFactors
#' @rdname normalize
#' @return Normalized data frame in slot normalized
#' @export
setGeneric(name = "NormalizespaceST",
           def = function(object,
                          method = "scran",
                          log2 = F,
                          clusters = NULL
           ) standardGeneric("NormalizespaceST"))
#' @rdname normalize
#' @export
setMethod(
  f = "NormalizespaceST",
  signature = "spaceST",
  definition = function(object,
                        method = "cp10k",
                        log2 = F,
                        clusters = NULL) {
    if (class(object) != "spaceST"){
      stop("Wrong input format.")
    }
    if (length(object@corrected) > 0){
      norm.data <- object@corrected
    } else {
      norm.data <- object@expr
    }
    if (log2 & method == "cp10k") {
      object@normalized <- log2(calc_cpm(norm.data) + 1)
    } else {
      if (method == "scran") {
        if (is.null(object@clusters))
          sce <- newSCESet(countData = data.frame(norm.data))
        if (!is.null(clusters)) {
          sce <- computeSumFactors(sce, cluster = clusters)
        } else {
          sce <- computeSumFactors(sce)
        }
        if (sum(sizeFactors(sce) == 0) > 0) {
          stop("Zero sum factors not allowed. Filter data before normalizetion with scran.")
        }
        sce <- scater::normalize(sce)
        object@normalized <- as.matrix(sce)
      } else if (method == "cp10k") {
        object@normalized <- calc_cpm(norm.data)
      }
    }
    return(object)
  }
)

#' Plot pca of spaceST data
#'
#' @description This function is used to plot the first two principal components of expression data
#' stored in a spaceST object. By default, both raw expression and batch corrected exåression
#' data will be plotted if both are present.
#' @param object spaceST object with expression data
#' @param ... arguments passed to pca_plot
#' @rdname pca
#' @return Plot of the first two principal components of expression data.
#' @export
setGeneric("pca.spaceST", function(object, ...) standardGeneric("pca.spaceST"))
#' @rdname pca
#' @export
setMethod("pca.spaceST", "spaceST", function(object, ...) {
  if (class(object) != "spaceST") {
    stop("Wrong input format")
  }
  if (length(object@corrected) == 0) {
    pca_plot(df1 = object@expr, df2 = NULL, samples = object@coordinates[, 1], ...)
  } else {
    pca_plot(object@expr, object@corrected, object@coordinates[, 1], ...)
  }
})

#' Plot unique genes per feature and transcripts per feature
#'
#' @description This funciton is used to plot the unique genes per feature and transcipts per feature
#' distributions as histograms. By default, the function will take the batch corrected dataset
#' if present.
#' @param object spaceST object with expression data.
#' @param dataset Character string specifying if the "corrected" or "raw" data should be used as input. Default
#' is set to "corrected".
#' @param separate Logical indicating whether replicates should be ploteed separately
#' @rdname QC
#' @return Histograms of unique genes per feature and transcripts per feature distributions.
#' @export
setGeneric("plot.QC.spaceST", function(object, dataset = "corrected", separate = F) standardGeneric("plot.QC.spaceST"))
#' @rdname QC
#' @docType methods
#' @export
setMethod("plot.QC.spaceST", "spaceST", function(object, dataset = "corrected", separate = F) {
  if (class(object) != "spaceST") {
    stop("Wrong input format")
  }
  if (length(object@corrected) == 0) {
    df1 <- ST_statistics(object@expr)
    ST_statistics.plot(df1, separate)
  }
  if (dataset == "corrected" & length(object@corrected) > 0) {
    df1 <- ST_statistics(object@corrected)
    ST_statistics.plot(df1, separate)
  }
  if (dataset == "raw") {
    df1 <- ST_statistics(object@expr)
    ST_statistics.plot(df1, separate)
  }
})
