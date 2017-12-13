#' xCellBoot class
#'
#' @slot rep.list  A list of pseudo-replicate bootstrap results
#' @slot obs.df Data.frame with observed xCell values
#' @slot mean.df Data.frame with mean bootstrap values
#' @slot bias.df Data.frame with bias values, i.e. difference between observed and mean values
#' @slot betadist p-values
#'
#' @name xCellBoot
#' @rdname xCellBoot
#' @aliases xCellBoot-class
#' @exportClass xCellBoot

xCellBoot <- setClass(
  "xCellBoot",
  slots = c(rep.list = "list",
            obs.df = "data.frame",
            mean.df = "data.frame",
            bias.df = "data.frame"
  )
)

#' Bootstrap xCell data.
#'
#' @description Used to bootstrap xCell analysis results.
#' @param data Input data.frame.
#' @param clusters Integer vector specifying cluster positions in data.frame.
#' @param replicates Integer specifying number of bootstrap replicates.
#' @return An xCellBoot object.
#' @importFrom xCell xCellAnalysis
#' @export
CreatexCellBoot <- function(
  data,
  clusters,
  replicates
  ) {
  if (!(class(data) %in% c("data.frame", "matrix"))) {
    stop("Wrong input format")
  }
  new(
    "xCellBoot"
    )
  object@rep.list <- xCell.bootstrap(
    data,
    clusters,
    replicates
  )
  object@obs.df <- as.data.frame(xCellAnalysis(t(rowsum(t(data), clusters))))
  celltypes <- rownames(object@obs.df)
  names(object@rep.list) <- celltypes
  rownames(object@obs.df) <- celltypes
  object@mean.df <- as.data.frame(do.call(rbind, lapply(X = object@rep.list,
                                                        FUN = function(x) (apply(x, 2, mean)))))
  rownames(object@mean.df) <- celltypes
  object@bias.df <- object@obs.df - object@mean.df
  rownames(object@bias.df) <- celltypes
  return(object)
}

#' Bootstrap function for xCell
#'
#' @description Calculate xCell scores for bootstrap replicates of clusters.
#' @param df Input gene expression data.frame.
#' @param clusters Integer vector specifying cluster identity for each feature.
#' @param replicates Integer specifying number of replicates to run.
#' @importFrom tcltk tkProgressBar
#' @importFrom xCell xCellAnalysis
#' @export
xCell.bootstrap <- function(df, clusters, replicates, progressbar = T) {
  if (progressbar) {
    rep_count <- 1
    # Suppress warnings globally
    options(warn = -1)
    pb <- tkProgressBar(min = 0,
                        max = replicates, width = 300)
  }
  rep.list <- list()
  for (n in 1:replicates) {
    sample.list <- list()
    clusters.new <- c()
    for (i in 1:max(unique(clusters))) {
      index <- (1:ncol(df))[clusters == i]
      sample.list[[i]] <- sample(index, length(index), replace = T)
      clusters.new <- c(clusters.new, rep(i, length(index)))
    }
    samples <- do.call(c, sample.list)
    sample.df <- df[, samples]
    clust.matrix <- t(rowsum(t(sample.df), clusters.new))
    xCell.res <- xCellAnalysis(clust.matrix)
    rep.list[[n]] <- xCell.res
    if (progressbar) {
      setTkProgressBar(pb, rep_count, label = paste(round(rep_count/replicates*100, 0), "% done"))
      rep_count <- rep_count + 1
    }
  }

  result.list <- list()
  for (i in 1:67) {
    result.list[[i]] <- do.call(rbind, lapply(rep.list, function(x) (x[i, ])))
  }
  close(pb)
  # Turn warning back on
  options(warn = 0)
  return(result.list)
}

setGeneric("celltype.bootstraps", function(object, celltype) standardGeneric("celltype.bootstraps"))

#' Extract celltype bootstrap replicates
#'
#' @description Extract the bootstrap results for a specific celltype.
#' @param celltype Character string specifying xCell signature celltype.
#' @return Data.frame with clusters as columns and replicates as rows.
#' @export
setMethod("celltype.bootstraps", "xCellBoot", function(object, celltype) {
  if (class(object) != "xCellBoot"){
    stop("Wrong input format.")
  }
  celltypes <- rownames(object@obs.df)
  index <- (1:67)[which(celltype %in% celltypes)]
  return(object@rep.list[[index]])
})

setGeneric("cluster.bootstraps", function(object, cluster, filter = F, probs = 0.1, threshold = 0.001) standardGeneric("cluster.bootstraps"))

#' Extract cluster bootstrap replicates
#'
#' @description Extract the bootstrap results for a specific cluster.
#' @param cluster Integer specifying cluster number.
#' #' @param filter Logical specifying wether or not filtering should be applied.
#' @param probs Chose quantile to use for filtering.
#' @param threshold Numeric value specifying threshold. If filter = TRUE, columns with probs_quantile > threshold will be kept.
#' @return Data.frame with celltypes as columns and replicates as rows.
#' @export
setMethod("cluster.bootstraps", "xCellBoot", function(object, cluster, filter = F, probs = 0.1, threshold = 0.001) {
  if (class(object) != "xCellBoot"){
    stop("Wrong input format.")
  }
  df <- do.call(cbind, lapply(object@rep.list, function(x) (x[, cluster])))
  if (filter) {
    indices.cells <- which(apply(df, 2, quantile, probs = probs) > threshold)
    df[, -indices.cells] <- NA
  }
  return(df)
})

setGeneric("obs.bootstraps", function(object, filter = F, probs = 0.1, threshold = 0.001) standardGeneric("obs.bootstraps"))

#' Extract obserevd bootstrap values
#'
#' @description Extract the observed xCell scores for each celltype and cluster from an xCellBoot object.
#' @param filter Logical specifying whether or not the data should be filtered from unreliable estimates.
#' @param probs Numerical value between 0-1 specifying the quantile used for filtering. (see threshold)
#' @param threshold Set threshold for xCell score for filtering. For a given celltype distribution of xCell scores,
#' the quantile value specified by the probs parameter will be used as a test to decide if the estimate should be kept or not.
#' The xCell score distributions might be heavily skewed with a large number of 0 values, which could arise due to the fact that
#' the celltype gene signature used to calculate those estimated are only present in a subset of the input features.
#' With the default filter settings, the 0.1 quantile is calculated for each distribution and if this value is larger
#' than the threshold of 0.001, the estimate remains unchanged. If the value is lower than the threshold, the
#' estimate is set to 0 instead. Distributions failing the test usually contain very low xCell scores but if the estimate
#' is too biased, this filtering serves as a means of omitting unreliable results.
#' @return Observed xCell scores for clusters.
#' @export
setMethod("obs.bootstraps", "xCellBoot", function(object, filter = F, probs = 0.1, threshold = 0.001) {
  if (class(object) != "xCellBoot"){
    stop("Wrong input format.")
  }
  df <- object@obs.df
  if (filter) {
    index.list <- list()
    for (i in 1:ncol(object@obs.df)) {
      index.list[[i]] <- which(apply(do.call(cbind, lapply(object@rep.list, function(x) (x[, i]))), 2, quantile, probs = probs) > threshold)
    }
    for (i in 1:ncol(df)) {
      df[-index.list[[i]], i] <- 0
    }
  }
  colnames(df) <- paste("c", colnames(df), sep = "")
  return(df)
})

setGeneric("mean.bootstraps", function(object) standardGeneric("mean.bootstraps"))

#' Extract mean bootstrap values
#'
#' @description Extract the mean xCell scores for each celltype and cluster from an xCellBoot object.
#' @return Data.frame with celltypes as columns and replicates as rows.
#' @export
setMethod("mean.bootstraps", "xCellBoot", function(object) {
  if (class(object) != "xCellBoot"){
    stop("Wrong input format.")
  }
  df <- object@mean.df
  colnames(df) <- paste("c", colnames(df), sep = "")
  return(df)
})

setGeneric("bias.bootstraps", function(object) standardGeneric("bias.bootstraps"))

#' Extract bias data
#'
#' @description Extract the mean xCell scores for each celltype and cluster from an xCellBoot object.
#' @return Data.frame with celltypes as columns and replicates as rows.
#' @export
setMethod("bias.bootstraps", "xCellBoot", function(object) {
  if (class(object) != "xCellBoot"){
    stop("Wrong input format.")
  }
  df <- object@bias.df
  colnames(df) <- paste("c", colnames(df), sep = "")
  return(df)
})

#' Violin plot of xCell scores
#'
#' @description Violin or boxplot of xCell scores.
#' @param df Input data.frame.
#' @param type Charcter string specifying plot type ('violin' or 'boxplot').
#' @importFrom reshape2 melt
#' @return Plot of xCell score distributions.
#' @export
violin.bootstraps <- function(df, type = "violin", cluster = NULL, group = NULL, samplename = NULL, ...) {
  if (!(type %in% c("violin", "boxplot"))) {
    stop("Wrong plot type.")
  }
  stopifnot(class(df) == "data.frame" | class(df) == "xCellBoot")
  if (!is.null(cluster) & class(df) != "xCellBoot") {
    stop("xCellBoot object required")
  }
  if (!is.null(cluster) & class(df) == "xCellBoot") {
    x <- cluster.bootstraps(object = df, cluster = cluster, ...)
    if (!is.null(group)) {
      x <- x[, which(colnames(x) %in% group)]
    }
    obs.data <- df@obs.df[, cluster][which(rownames(df@obs.df) %in% colnames(x))]
    df <- x
  }
  dfm <- as.data.frame(melt(df))
  colnames(dfm) <- c("id", "names", "val")
  p <- ggplot2::ggplot(dfm, ggplot2::aes(names, val, fill = names))
  if (type == "boxplot") {
    p <- p + ggplot2::stat_boxplot(geom ='errorbar', width = 0.5) +
      ggplot2::geom_boxplot()
  } else if (type == "violin") {
    p <- p + ggplot2::geom_violin(scale = "width") +
      ggplot2::geom_boxplot(width = .1)
  }
  if (!is.null(cluster)) {
    obs.data <- data.frame(names = colnames(x), obs.data)
    p <- p + ggplot2::geom_point(data = obs.data, ggplot2::aes(names, obs.data), stat = "identity", color = "red")
    if (!is.null(samplename)) {
      p <- p + ggplot2::labs(title = samplename)
    }
  }
  p <- p + ggplot2::labs(x = "celltype", y = "xCell score") +
    ggplot2::guides(fill = "none") +
    ggplot2::theme(panel.background = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(colour = "black", fill = NA, size = 0.2),
                                  axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))
  plot(p)
}
