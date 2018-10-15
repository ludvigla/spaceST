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
    status.expr = "character",
    DE.list = "ANY"
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
  cat("  Filtered out spots with less than ", object@filter.settings[[1]], " unique genes \n", sep = "")
  cat("  Filtered out genes with expression values =< ", object@filter.settings[[2]], " in =< ", object@filter.settings[[3]],
      " features\n", sep = "")
  cat("  Filtered out genes matching ", object@filter.settings[[4]], " \n\n", sep = "")
  cat("expr status:", object@status.expr)
})


#' Method batch.correct for spaceST class.
#'
#' Function used to run a batch correction between samples.
#' @param object Object of class spaceST.
#' @seealso \link[CountClust]{BatchCorrectedCounts}
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
#' @name NormalizespaceST
#' @rdname spaceST-methods
#' @param object Object of class spaceST.
#' @param method Select normalization method, default "cp10k" [options: "cp10k", "scran"].
#' @param log2 Logical specifying whether or not data should be log2-transformed.
#' @param pcount Pseudocount value for log2 tranformation (default = 1).
#' @param clusters Integer vector specifying clusters used for normalization with scran on large gene expression matrices.
#' By default, a clustering step is performed using \link[scran]{quickCluster} if there are more than 500 array spots in the
#' expression matrix.
#' @seealso \link[scater]{normalize} and \link[scran]{computeSumFactors}
#' @importFrom scran quickCluster computeSumFactors
#' @importFrom SingleCellExperiment logcounts sizeFactors
#' @return Normalized expression data.
#' @exportMethod NormalizespaceST
setGeneric(name = "NormalizespaceST",
           def = function(object,
                          method = "scran",
                          log2 = F,
                          pcount = 1,
                          clusters = NULL
           ) standardGeneric("NormalizespaceST"))
#' @rdname spaceST-methods
#' @aliases NormalizespaceST,spaceST-methods
setMethod(
  f = "NormalizespaceST",
  signature = "spaceST",
  definition = function(object,
                        method = "scran",
                        log2 = F,
                        pcount = 1,
                        clusters = NULL) {
    if (class(object) != "spaceST"){
      stop("Wrong input format.")
    }
    expr.data <- as.matrix(object@expr)
    if (log2 & method == "cp10k") {
      object@norm.data <- log2(calc_cpm(expr.data) + 1)
    } else {
      if (method == "scran") {
        sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = expr.data, logcounts = log2(expr.data + pcount)))
        if (ncol(object@expr) > 700 | !is.null(clusters)) {
          if (!is.null(clusters)) {
            sce <- scran::computeSumFactors(sce, cluster = clusters)
            print("Using provided clusters ...")
          } else {
            q.clust <- quickCluster(as.matrix(expr.data))
            sce <- scran::computeSumFactors(sce, cluster = q.clust)
          }

          object@meta.data$size_factor <- SingleCellExperiment::sizeFactors(sce)
        } else {
          sce <- computeSumFactors(sce)
          object@meta.data$size_factor <- SingleCellExperiment::sizeFactors(sce)
        }
        if (sum(sizeFactors(sce) == 0) > 0) {
          stop("Zero sum factors not allowed. Filter data before normalization with scran.")
        }
        sce <- scater::normalize(sce)
        object@norm.data <- as(SingleCellExperiment::logcounts(sce), "dgCMatrix")
      } else if (method == "cp10k") {
        object@norm.data <- as(calc_cp10k(as.matrix(expr.data)), "dgCMatrix")
      }
    }
    return(object)
  }
)


#' Method plotPCA.
#'
#' Plot principal components.
#' @name plotPCA
#' @rdname spaceST-methods
#' @param object Object of class spaceST.
#' @param components Integer or character vector of length 2 specifying two PC components to compare.
#' @param ... Parameters passed to geom_point().
#' @return PCA plot.
#' @importFrom ggplot2 ggplot aes_string geom_point theme_classic scale_color_brewer labs
#' @exportMethod plotPCA
setGeneric(name = "plotPCA",
           def = function(object,
                          components = c(1, 2),
                          ...
           ) standardGeneric("plotPCA"))
#' @rdname spaceST-methods
#' @aliases plotPCA,spaceST-methods
setMethod(
  f = "plotPCA",
  signature = "spaceST",
  definition = function(object,
                        components = c(1, 2),
                        ...
                        ) {

  if (length(object@reducedDims) == 0) {
    message("reducedDims not present in spaceST object. Running PC analysis using normalized data ...")
    stopifnot(length(object@norm.data) == 0)
    object <- spPCA(object)
  }
  pcs <- as.data.frame(object@reducedDims$x)
  pov <- summary(object@reducedDims)$importance
  pcs$replicate <- object@coordinates$replicate
  if (class(components) %in% c("integer", "numeric")) {
    components <- colnames(pcs)[components]
  }
  ggplot(pcs, aes_string(components[1], components[2], color = "replicate")) +
    geom_point(...) +
    theme_classic() +
    scale_color_brewer(palette = "Dark2") +
    labs(x = paste(components[1], round(pov[2, components[1]]*100, digits = 2), "%"), y = paste(components[2], round(pov[2, components[2]]*100, digits = 2), "%"))
})



#' Method plotTSNE.
#'
#' Plot t-SNE.
#' @name plotTSNE
#' @rdname spaceST-methods
#' @param object Object of class spaceST.
#' @param group_by character; Grouping variable for color.
#' @param cols character; Specify colors for grouping variable.
#' @param ... Parameters passed to geom_point().
#' @return t-SNE plot.
#' @importFrom ggplot2 ggplot aes_string geom_point theme_classic scale_color_brewer labs
#' @exportMethod plotTSNE
setGeneric(name = "plotTSNE",
           def = function(object,
                          group_by = NULL,
                          cols = NULL,
                          ...
           ) standardGeneric("plotTSNE"))
#' @rdname spaceST-methods
#' @aliases plotPCA,spaceST-methods
setMethod(
  f = "plotTSNE",
  signature = "spaceST",
  definition = function(object,
                        group_by = NULL,
                        cols = NULL,
                        ...
  ) {

    if (length(object@tsne) == 0) {
      stop("tsne not present in spaceST object ...")
    }

    gg.df <- as.data.frame(object@tsne[, 1:2])
    colnames(gg.df) <- c("x", "y")

    if (is.null(cols)) {
      cols <- c("royalblue3", "mediumpurple4", "navajowhite2", "chocolate", "firebrick", "yellow2"," aquamarine", "orange1", "olivedrab2", "darkgreen", "pink", "black", "navy", "khaki3", "lightsteelblue1")
    }

    if (is.null(group_by)) {
      p <- ggplot(gg.df, aes(x = x, y = y))
    } else if (group_by == "clusters") {
      clusters <- object@meta.data$clusters
      p <- ggplot(gg.df, aes(x = x, y = y, colour = factor(clusters))) +
        scale_color_manual(values = cols) +
        labs(colour = "cluster", title = paste("t-SNE", sep = ""))
    } else if (group_by == "clusters.snn") {
      clusters <- object@meta.data$clusters.snn
      p <- ggplot(gg.df, aes(x = x, y = y, colour = factor(clusters))) +
        scale_color_manual(values = cols) +
        labs(colour = "cluster", title = paste("t-SNE", sep = ""))
    } else if (group_by == "sample") {
      sample <- object@coordinates$replicate
      p <- ggplot(gg.df, aes(x = x, y = y, colour = factor(sample))) +
        scale_color_manual(values = cols) +
        labs(colour = "sample", title = paste("t-SNE", sep = ""))
    }

    p <- p + geom_point(...) +
      theme_classic() +
      labs(x = "t-SNE dim 1", y = "t-SNE dim 2")
    plot(p)
    })




#' Plot unique genes per feature and transcripts per feature.
#'
#' This function is used to plot the unique genes per feature and transcipts per feature
#' distributions as histograms. By default, the function will take the batch corrected dataset
#' if present.
#' @param object Object of class spaceST.
#' @param separate Logical indicating whether replicates should be ploteed separately.
#' @param bins Select number of bins.
#' @param ... Parameters passed to geom_histogram.
#' @return Histograms of unique genes per feature and transcripts per feature distributions.
#' @export
setGeneric("plot_QC_spaceST", function(object, separate = F, bins = 20, type = "norm.data", ...) standardGeneric("plot_QC_spaceST"))
#' @name plot_QC_spaceST
#' @aliases plot_QC_spaceST,spaceST-methods
#' @rdname spaceST-methods
#' @export
setMethod(f = "plot_QC_spaceST", signature = "spaceST", definition = function(object, separate = F, bins = 20, type = "norm.data", ...) {
  if (class(object) != "spaceST") {
    stop("Wrong input format")
  }
  samples <- object@coordinates$replicate
  if (length(object@norm.data) > 0 & type == "norm.data") {
    df <- cbind(ST_statistics(as.matrix(object@norm.data)), samples = samples)
    ST_statistics_plot(df, separate, bins = bins, ...)
  } else if (type == "expr") {
    df <- cbind(ST_statistics(as.matrix(object@expr)), samples = samples)
    ST_statistics_plot(df, separate, bins = bins, ...)
  }
})


#' Method spPCA for spaceST class.
#'
#' Run PC analysis on spaceST object.
#' @param object Object of class spaceST.
#' @param ntop Number of variable genes to select for PC analysis.
#' @param exprs_values Select target object to use as input for PC analysis [options: "norm.data", "expr"]
#' @param log2.transform Log2-transform data before running pca.
#' @param ... Parameters passed to \link[stats]{prcomp}.
#' @seealso \link[stats]{prcomp}
#' @export
setGeneric("spPCA", function(object,
                             ntop = 500,
                             exprs_values = "norm.data",
                             log2.transform = TRUE,
                             ...) standardGeneric("spPCA"))
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
  log2.transform = TRUE,
  ...
) {
  if (exprs_values == "norm.data") {
    input.data <- as.matrix(object@norm.data)
  } else if (exprs_values == "exprs") {
    input.data <- as.matrix(object@expr)
  }
  if (log2.transform) {
    input.data <- log2(input.data)
  }
  high.genes <- names(sort(apply(input.data, 1, var), decreasing = T)[1:ntop])
  input.data <- input.data[high.genes, ]
  pca <- prcomp(t(input.data), scale. = T, center = T, ...)

  object@reducedDims <- pca
  return(object)
})



#' Method spots_under_tissue for spaceST class.
#'
#' Subset spaceST object using spot selection tables from the ST spot detector.
#' @param object Object of class spaceST.
#' @param selection.files List of paths to selection files. The order of selections has to match the order of
#' sample matrices used to initiate the spaceST object.
#' @param delimiter Delimiter for expr headers.
#' @param keep.filter Set to FALSE if you don't want to apply the same filtering settings that were used when the
#' spaceST object was initiated.
#' @seealso \link[CreatespaceSTobject]{spaceST}
#' @export
setGeneric("spots_under_tissue", function(object, selection.files, delimiter = "_", keep.filter = T) standardGeneric("spots_under_tissue"))
#' @rdname spaceST-methods
#' @name spots.under.tissue
#' @aliases spots.under.tissue,spaceST-methods
#' @export
setMethod(
  "spots_under_tissue",
  signature = "spaceST",
  definition = function(
    object,
    selection.files,
    delimiter = "_",
    keep.filter = T
) {
  stopifnot(class(object) == "spaceST")
  reps <- unique(object@coordinates$replicate)
  all.coords <- object@coordinates
  expr <- as.matrix(object@expr)
  selection.list <- lapply(1:length(reps), function(i) {
    # Select replicate subset
    expr.subset <- expr[, all.coords$replicate == reps[i]]
    # Select replicate coordinates
    coords <- all.coords[all.coords$replicate == reps[i], 2:3]
    # Create character vector of selected coordinates
    spots <- paste(round(coords[, 1]), round(coords[, 2]), sep = "x")
    # Read alignment table
    alignment <- read.table(selection.files[[i]], header = T)
    if (ncol(alignment) == 7) {
      alignment <- alignment[alignment[, 7] == 1, ]
    }
    # Create character vector of spots under tissue
    alignment.spots <- paste(alignment$x, alignment$y, sep = "x")
    # Sanity check
    stopifnot(sum(alignment.spots %in% spots) == sum(spots %in% alignment.spots))
    # Intersecting spots
    intersecting.spots <- which(alignment.spots %in% spots)

    # Select iondices for spots under tissue
    spot_indices <- which(spots %in% alignment.spots)
    # Create subset of input matrix with spots under tissue
    subset_expr <- expr.subset[, spot_indices]
    # Subset alignment table with intersecting spots
    alignment <- alignment[intersecting.spots, ]
    # Create dictionary of new coordinates
    new.spots <- paste(reps[i], round(alignment$new_x, digits = 2), round(alignment$new_y, digits = 2), sep = delimiter)
    names(new.spots) <- paste(reps[i], round(alignment$x, digits = 2), round(alignment$y, digits = 2), sep = delimiter)
    colnames(subset_expr) <- new.spots[colnames(subset_expr)]
    return(subset_expr)
  })
  expr <- do.call(cbind, selection.list)
  if (keep.filter) {
    new.spST <- CreatespaceSTobject(expr, delimiter = delimiter,
                                    unique.genes = object@filter.settings$unique.genes,
                                    min.features = object@filter.settings$min.features,
                                    min.exp = object@filter.settings$min.exp,
                                    filter.genes = object@filter.settings$filter.genes
                                    )
  } else {
    new.spST <- CreatespaceSTobject(expr, delimiter = delimiter)
  }
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
#' @param force.filter Force filter on a spaceST object.
#' @param delimiter Character used to delimit coordinates in headers.
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
    object <- CreatespaceSTobject(raw.data = as.matrix(object@expr), unique.genes, min.exp, min.features, filter.genes, delimiter)
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
    object@coordinates <- get_coordinates(expr)
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
#' @param object Object of class spaceST.
#' @exportMethod counts
setGeneric("counts", function(object) standardGeneric("counts"))

#' @name counts
#' @rdname spaceST-methods
#' @param object Object of class spaceST.
#' @aliases counts,spaceST-methods
setMethod("counts", "spaceST", function(object) as.matrix(object@expr))

#' Method normcounts.
#' @name normcounts
#' @param object Object of class spaceST.
#' @exportMethod normcounts
setGeneric("normcounts", function(object) standardGeneric("normcounts"))

#' @name normcounts
#' @rdname spaceST-methods
#' @param object Object of class spaceST.
#' @aliases normcounts,spaceST-methods
setMethod("normcounts", "spaceST", function(object) as.matrix(object@norm.data))

#' Method topic.clusters.
#' @name topic.clusters
#' @param @param object Object of class spaceST.
#' @exportMethod topic.clusters
setGeneric("topic.clusters", function(object) standardGeneric("topic.clusters"))

#' @name topic.clusters
#' @rdname spaceST-methods
#' @param object Object of class spaceST.
#' @aliases topic.clusters,spaceST-methods
setMethod("topic.clusters", "spaceST", function(object) object@meta.data$clusters)

#' Method topic.clusters.
#' @name getLDA
#' @param object Object of class spaceST.
#' @exportMethod getLDA
setGeneric("getLDA", function(object) standardGeneric("getLDA"))

#' @name getLDA
#' @rdname spaceST-methods
#' @param object Object of class spaceST.
#' @aliases getLDA,spaceST-methods
setMethod("getLDA", "spaceST", function(object) object@lda.results)
