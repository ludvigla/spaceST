#' t-SNE of spaceST data
#'
#' @description Run t-SNE on spaceST object. Results are saved in the tsne slot.
#' @param object A spaceST object
#' @param dims integer; Output dimensionality (default: 2)
#' @param initial.dims integer; the number of dimensions that should be retained in the initial PCA step (default: 50).
#' Only applicable if pca = TRUE.
#' @param theta numeric; Speed/accuracy trade-off (increase for less accuracy), set to 0.0 for exact TSNE (default: 0.5).
#' @param check_duplicates logical; Checks whether duplicates are present. It is best to make sure there are no duplicates present
#' and set this option to FALSE, especially for large datasets (default: TRUE).
#' @param pca logical; Whether an initial PCA step should be performed (default: FALSE.
#' @param perplexity numeric; Perplexity parameter.
#' @param max_iter integer; Number of iterations (default: 1000).
#' @param verbose logical; Whether progress updates should be printed (default: FALSE).
#' @param select.dims integer; Select specific dimensions of input data, i.e. columns (default: NULL).
#' @param seed Set seed for reproducibility.
#' @param plot.tSNE logical specifying whether or not t-SNE results should be plotted.
#' @param use.dims specify what dimensions should be used for t-SNE (currently only "topics" available).
#' @param clusters integer vector used to color features in t-SNE plot.
#' @param cols character; Specify colors for grouping variable.
#' @param ... Parameters passed to Rtsne.
#' @seealso \link[Rtsne]{Rtsne}
#' @importFrom Rtsne Rtsne
#' @importFrom graphics plot
#' @return Matrix with t-SNE results
#' @rdname RunTSNEspaceST
#' @export
RunTSNEspaceST <- function(
  object,
  dims = 2,
  initial.dims = 50,
  theta = 0.0,
  check_duplicates = FALSE,
  pca = FALSE,
  perplexity = 15,
  max_iter = 1000,
  verbose = FALSE,
  seed = 0,
  plot.tSNE = FALSE,
  use.dims = "topics",
  select.dims = NULL,
  group_by = NULL,
  cols = NULL,
  ...
) {
  if (use.dims == "topics") {
    stopifnot(length(object@lda.results) > 0)
    df <- object@lda.results$omega
  } else if (use.dims == "pca") {
    stopifnot(length(object@reducedDims) > 0)
    df <- object@reducedDims$x
    if (is.null(select.dims)) {
      df <- df[, select.dims]
    }
  } else if (use.dims == "expr") {
    stopifnot(length(object@expr) > 0)
    df <- as.matrix(object@expr)
    if (is.null(select.dims)) {
      vars <- sort(apply(df, 1, var), decreasing = T)
      if (length(vars) > 500) {
        vars <- vars[1:500]
      }
      df <- t(df[names(vars), ])
    } else {
      df <- t(df[select.dims, ])
    }
  } else if (use.dims == "norm.data") {
    stopifnot(length(object@norm.data) > 0)
    df <- as.matrix(object@norm.data)
    if (is.null(select.dims)) {
      vars <- sort(apply(df, 1, var), decreasing = T)
      if (length(vars) > 500) {
        vars <- vars[1:500]
      }
      df <- t(df[names(vars), ])
    } else {
      df <- t(df[select.dims, ])
    }
  }
  # TODO: addmethods
  set.seed(seed)
  out.tsne = Rtsne(as.matrix(df),
                   dims = dims,
                   initial_dims = initial.dims,
                   theta = theta,
                   check_duplicates = check_duplicates,
                   pca = pca,
                   perplexity = perplexity,
                   max_iter = max_iter,
                   verbose = verbose,
                   ...)
  if (plot.tSNE) {
    if (!is.null(clusters)) {
      gg.df <- data.frame(x = out.tsne$Y[,1], y = out.tsne$Y[,2], clusters = clusters)
    } else {
      gg.df <- data.frame(x = out.tsne$Y[,1], y = out.tsne$Y[,2])
    }

    if (is.null(cols)) {
      cols <- c("royalblue3", "mediumpurple4", "navajowhite2", "chocolate", "firebrick", "yellow2"," aquamarine", "orange1", "olivedrab2", "darkgreen", "pink", "black", "navy", "khaki3", "lightsteelblue1")
    }
    if (is.null(group_by)) {
      p <- ggplot(gg.df, aes(x = x, y = y))
    } else if (group_by == "clusters") {
      clusters <- object@meta.data$clusters
      p <- ggplot(gg.df, aes(x = x, y = y, colour = factor(clusters))) +
        scale_color_manual(values = cols) +
        labs(colour = "cluster", title = paste("t-SNE of ", use.dims, " results", sep = ""))
    } else if (group_by == "clusters.snn") {
      clusters <- object@meta.data$clusters.snn
      p <- ggplot(gg.df, aes(x = x, y = y, colour = factor(clusters))) +
        scale_color_manual(values = cols) +
        labs(colour = "cluster", title = paste("t-SNE of ", use.dims, " results", sep = ""))
    } else if (group_by == "sample") {
      sample <- object@coordinates$replicate
      p <- ggplot(gg.df, aes(x = x, y = y, colour = factor(sample))) +
        scale_color_manual(values = cols) +
        labs(colour = "sample", title = paste("t-SNE of ", use.dims, " results", sep = ""))
    }

    p <- p + geom_point() +
      theme(panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            panel.border = element_rect(colour = "black", fill = NA, size = 1))
    plot(p)
  }
  object@tsne <- out.tsne$Y
  return(object)
}
