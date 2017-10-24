#' t-SNE of spaceST data
#'
#' @description Run t-SNE on spaceST object. Results are saved in the tsne slot.
#' @param object A spaceST object
#' @param seed Set seed for reproducibility
#' @param plot.tSNE Logical specifying whether or not t-SNE results should be plotted
#' @param use.dims Specify what dimensions should be used for t-SNE (currently only "topics" available)
#' @param clusters Integer vector used to color features in t-SNE plot
#' @importFrom Rtsne Rtsne
#' @return Matrix with t-SNE results
#' @rdname RunTSNEspaceST
#' @export
RunTSNEspaceST <- function(
  object,
  dims = 2,
  initial.dims = 50,
  theta = 0.0,
  check_duplicates = FALSE,
  pca = TRUE,
  perplexity = 15,
  max_iter = 1000,
  verbose = FALSE,
  seed = 0,
  plot.tSNE = FALSE,
  use.dims = "topics",
  clusters = NULL
) {
  if (use.dims == "topics") {
    stopifnot(length(object@topics) > 0)
    df <- object@topics
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
                   verbose = verbose)
  if (plot.tSNE) {
    if (!is.null(clusters)) {
      gg.df <- data.frame(x = out.tsne$Y[,1], y = out.tsne$Y[,2], clusters = clusters)
    } else {
      gg.df <- data.frame(x = out.tsne$Y[,1], y = out.tsne$Y[,2])
    }
    cols <- c("royalblue3", "mediumpurple4", "navajowhite2", "chocolate", "firebrick", "yellow2"," aquamarine", "orange1", "olivedrab2", "darkgreen", "pink", "black", "navy", "khaki3", "lightsteelblue1")
    if (!is.null(clusters)) {
      p <- ggplot(gg.df, aes(x = x, y = y, colour = factor(clusters))) +
        scale_color_manual(values = cols) +
        labs(colour = "cluster", title = paste("t-SNE of ", dims.use, " results", sep = ""))
    } else {
      p <- ggplot(gg.df, aes(x = x, y = y))
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
