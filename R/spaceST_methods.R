#' space-ST methods.
#' @import methods
#' @param object Object of class spaceST.
#' @name spaceST-methods
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


#' Method show for spaceST class.
#'
#' @param object Object of class spaceST.
#' @aliases show,spaceST-methods
#' @name show
#' @rdname spaseST-methods
setMethod("show", signature = "spaceST", definition = function(object) {
  cat("An object of class ", class(object), "\n\n", sep = "")
  cat("Filter settings:\n")
  cat("  Filtered features with less than ", object@filter.settings[[1]], " unique genes \n", sep = "")
  cat("  Filtered genes with expression values =< ", object@filter.settings[[2]], " in =< ", object@filter.settings[[3]],
      " features\n\n", sep = "")
  cat("expr status:", object@status.expr)
})


#' Method batch.correct for spaceST class.
#'
#' Function used to run a batch correction between samples.
#' @param object Object of class spaceST.
#' @seealso \link[countClust]{BatchCorrectedCounts}
#' @importFrom CountClust BatchCorrectedCounts
#' @return Batch corrected data frame in slot "expr".
#' @name batch.correct
#' @export
setGeneric(
  name = "batch.correct",
  def = function(object) standardGeneric("batch.correct")
)
#' @name batch.correct
#' @aliases batch.correct,spaceST-methods
#' @rdname spaceST-methods
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


#' Method NormalizespaceST ST data.
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
#' @return Normalized expression data.
#' @export
setGeneric(name = "NormalizespaceST",
           def = function(object,
                          method = "scran",
                          log2 = F,
                          clusters = NULL, ...
           ) standardGeneric("NormalizespaceST"))
#' @name NormalizespaceST
#' @aliases NormalizespaceST,spaceST-methods
#' @rdname spaceST-methods
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
#' This function is used to plot the unique genes per feature and transcipts per feature
#' distributions as histograms. By default, the function will take the batch corrected dataset
#' if present.
#' @param object Object of class spaceST.
#' @param separate Logical indicating whether replicates should be ploteed separately.
#' @return Histograms of unique genes per feature and transcripts per feature distributions.
#' @export
setGeneric("plot.QC.spaceST", function(object, separate = F) standardGeneric("plot.QC.spaceST"))
#' @name plot.QC.spaceST
#' @aliases plot.QC.spaceST,spaceST-methods
#' @rdname spaceST-methods
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


#' Method spPCA for spaceST class.
#'
#' Run PC analysis on spaceST object.
#' @param object Object of class spaceST.
#' @param ntop Number of variable genes to select for PC analysis.
#' @param ncomponents Number of components that should be returned.
#' @param exprs_values Select target object to use as input for PC analysis [options: "norm.data", "expr"]
#' @param ... Parameters passed to \link[stats]{prcomp}.
#' @seealso \link[stats]{prcomp}
#' @export
setGeneric("spPCA", function(object, ntop = 500, ncomponents = 2, exprs_values = "norm.data", ...) standardGeneric("spPCA"))
#' @export
#' @name spPCA
#' @aliases spPCA,spaceST-methods
#' @rdname spaceST-methods
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


#' Method spots.under.tissue for spaceST class.
#'
#' Subset spaceST object using spot selection tables from the ST spot detector.
#' @param object Object of class spaceST.
#' @param selection.files List of paths to selection files. The order of selections has to match the order of
#' sample matrices used to initiate the spaceST object.
#' @param delimiter Delimiter for expr headers.
#' @seealso \link[CreatespaceSTobject]{spaceST}
#' @export
setGeneric("spots.under.tissue", function(object, selection.files, delimiter = "_") standardGeneric("spots.under.tissue"))
#' @rdname spaceST-methods
#' @name spots.under.tissue
#' @aliases spots.under.tissue,spaceST-methods
#' @export
setMethod(
  "spots.under.tissue",
  signature = "spaceST",
  definition = function(
    object,
    selection.files,
    delimiter = "_"
) {
  stopifnot(class(object) == "spaceST")
  reps <- unique(object@coordinates$replicate)
  all.coords <- object@coordinates
  expr <- as.matrix(object@expr)
  selection.list <- lapply(1:length(reps), function(i) {
    coords <- all.coords[all.coords$replicate == reps[i], 2:3]
    spots <- paste(round(coords[, 1]), round(coords[, 2]))
    alignment <- read.table(selection.files[[i]], header = T)
    if (ncol(alignment) == 7) {
      alignment <- subset(alignment, selected == 1)
    }
    alignment.spots <- paste(alignment$x, alignment$y)
    stopifnot(sum(alignment.spots %in% spots) == sum(spots %in% alignment.spots))
    intersecting.spots <- which(alignment.spots %in% spots)
    spot_indices <- which(spots %in% alignment.spots)
    subset_expr <- expr[, spot_indices]
    alignment <- alignment[intersecting.spots, ]
    colnames(subset_expr) <- paste(i, round(alignment$new_x, digits = 2), round(alignment$new_y, digits = 2), sep = delimiter)
    return(subset_expr)
  })
  expr <- do.call(cbind, selection.list)
  new.spST <- CreatespaceSTobject(expr, delimiter = delimiter)
  return(new.spST)
})


#' Method filter for spaceST class.
#'
#' Apply filter to spaceST data.
#' @param object Object of class spaceST.
#' @param unique.genes The lowest number of unique genes allowed in a feature (spot).
#' @param min.exp Integer value specifying the lowest expression value allowed at least min.features number of features.
#' @param min.features Integer value specifying the lowest number of features with at least min.exp.
#' @param filter.genes A character vector specifying genes that should be filtered from the expression data.
#' @param delimiter Delimiter specifying header format.
#' @export
setGeneric("filter", function(object, unique.genes = 300, min.exp = 2, min.features = 15, filter.genes = NULL, force.filter = FALSE, delimiter = "_") standardGeneric("filter"))
#' @name filter
#' @aliases filter,spaceST-methods
#' @rdname spaceST-methods
#' @export
setMethod(
  "filter",
  signature = "spaceST",
  definition = function(
    object,
    unique.genes = 300,
    min.exp = 2,
    min.features = 15,
    filter.genes = NULL,
    force.filter = FALSE,
    delimiter = "_"
) {
  check <- any(length(object@norm.data) > 0,
               length(object@lda.results) > 0,
               length(object@reducedDims) > 0,
               length(object@tsne) > 0)
  if (check & !force.filter) {
    stop("Set force.filter = TRUE to initiate new spaceST object with supplied filtering settings. All data excep 'expr' will be lost.")
  } else if (check & force.filter) {
    object <- CreatespaceSTobject(as.matrix(object@expr), unique.genes, min.exp, min.features, filter.genes, delimiter)
  } else {
    expr <- object@expr
    if (!is.null(filter.genes)){
      removed.genes <- grep(rownames(expr), pattern = filter.genes, perl = TRUE)
      if(length(removed.genes) > 0) {
        expr = expr[-removed.genes, ]
      }
    }
    # Filter out low quality genes
    expr = expr[Matrix::rowSums(expr >= min.exp) >= min.features, ]

    # Filter out low quality spots
    indices <- which(apply(expr, 2, function(x) sum(x > 0)) < unique.genes)
    if (length(indices) > 0){
      expr <- expr[, -indices]
    }

    filter.settings = list(
      unique.genes = unique.genes,
      min.exp = min.exp,
      min.features = min.features,
      filter.genes = filter.genes)

    object@expr <- expr
    object@filter.settings <- filter.settings
    return(object)
  }
})

#' Method dim
#' @name dim
#' @rdname spaceST-methods
#' @aliases dim,spaceST-methods
setMethod("dim", "spaceST", function(x) dim(x@expr))

#' Method counts.
#' @name counts
#' @exportMethod counts
setGeneric("counts", function(object) standardGeneric("counts"))

#' @name counts
#' @rdname spaceST-methods
#' @aliases counts,spaceST-methods
setMethod("counts", "spaceST", function(object) as.matrix(object@expr))

#' Method normcounts.
#' @name normcounts
#' @exportMethod normcounts
setGeneric("normcounts", function(object) standardGeneric("normcounts"))

#' @name normcounts
#' @rdname spaceST-methods
#' @aliases normcounts,spaceST-methods
setMethod("normcounts", "spaceST", function(object) as.matrix(object@norm.data))

#' Method topic.clusters.
#' @name topic.clusters
#' @exportMethod topic.clusters
setGeneric("topic.clusters", function(object) standardGeneric("topic.clusters"))

#' @name topic.clusters
#' @rdname spaceST-methods
#' @aliases topic.clusters,spaceST-methods
setMethod("topic.clusters", "spaceST", function(object) as.matrix(object@meta.data$clusters))
