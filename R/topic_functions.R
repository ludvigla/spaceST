#' Helper function to compute lda using cellTree package.
#'
#' This helper function is used to compute topics for an expression dataset using lda modelling from cellTree package.
#' @param object Object of class spaceST.
#' @param min.topics,max.topics Integer values specifying min and max number of topics in the "maptpx" model.
#' @param max.topic Maximum number of topics.
#' @param num.genes Number of genes to use for lda modelling.
#' @param datatype Character string specifying if the "expr" or "norm.data" ata should be used as input. Default
#' is set to "norm.data".
#' @param method Method passed to CellTree function.
#' @param sd.filter Standard deviation threshold.
#' @param log.scale Convert to log scale.
#' @seealso \link[cellTree]{compute.lda}
#' @references \url{https://bioconductor.org/packages/release/bioc/html/cellTree.html}
#' @importFrom cellTree compute.lda
#' @return Topic matrix
#' @export
topic_compute <- function(
  object,
  min.topic = 2,
  max.topic = 15,
  num.genes = NULL,
  datatype = "norm.data",
  method = "maptpx",
  sd.filter = F,
  log.scale = F,
  force.recalc = F,
  minClusterSize = 30
) {
  UseMethod("topic_compute")
}
#' @export
topic_compute.default <- function(
  object,
  min.topic = 2,
  max.topic = 15,
  num.genes = NULL,
  method = "maptpx",
  sd.filter = F,
  log.scale = F,
  force.recalc = F,
  minClusterSize = 30
){
  ave.exp = rowMeans(object)
  sort.ave.exp = sort(ave.exp, decreasing = T)
  if (!is.null(num.genes)) {
    high.ave.exp = sort.ave.exp[c(1:num.genes)]
    high.genes = names(high.ave.exp)
  } else {
    high.genes = names(sort.ave.exp)
  }
  GoM.sample.df = object[high.genes,]
  K = c(min.topic:max.topic)
  lda.results = compute.lda(GoM.sample.df,
                            k.topics = K,
                            method = method,
                            sd.filter = sd.filter,
                            log.scale = log.scale)
  return(lda.results)
}
#' @export
topic_compute.spaceST <- function(
  object,
  min.topic = 2,
  max.topic = 15,
  num.genes = NULL,
  datatype = "norm.data",
  method = "maptpx",
  sd.filter = FALSE,
  log.scale = FALSE,
  force.recalc = FALSE,
  minClusterSize = 30
){
  if (!length(object@lda.results) > 0 | force.recalc) {
    if (force.recalc) {
      warning("Overwriting LDA results.")
    }
    object@meta.data <- list()
  } else if (!length(object@lda.results) > 0 | !force.recalc) {
    stop("LDA method has already been computed for this object. Set force.calc = TRUE if you want to overwrite the results")
  }
  if (datatype == "norm.data") {
    df <- as.matrix(object@norm.data)
  } else if (datatype == "expr") {
    df <- as.matrix(object@expr)
  }
  lda.results <- topic_compute.default(
    df,
    min.topic = 2,
    max.topic = 15,
    num.genes = NULL,
    method = "maptpx",
    sd.filter = F,
    log.scale = F
  )
  object@lda.results <- lda.results
  object@meta.data$clusters <- topic_clusters(omega = lda.results$omega, minClusterSize = minClusterSize)
  object@meta.data$minClusterSize <- minClusterSize
  return(object)
}


#' Calculate clusters based on topics.
#'
#' This function is used to cluster features based on a topics matrix.
#' @param omega Topic data.frame.
#' @param method.dist Set distance method.
#' @param method.tree Set clustering method.
#' @param minClusterSize Integer value specifying the minimum cluster size allowed.
#' @return Integer vector specifying cluster identity of each feature.
#' @export
topic_clusters <- function(omega, method.dist = "euclidean", method.tree = "ward.D2", minClusterSize = 30){
  my.dist = dist(omega, method = method.dist)
  my.tree = hclust(my.dist, method = method.tree)
  clusters = unname(dynamicTreeCut::cutreeDynamic(my.tree, distM = as.matrix(my.dist), verbose = 0, minClusterSize = minClusterSize))
  return(clusters)
}


#' Plot heatmap of lda results.
#'
#' This function is used to plot a heatmap of lda results.
#' @param df Topic data.frame (omega) or "topics" class object.
#' @param method.dist Set distance method.
#' @param method.tree Set clustering method.
#' @param minClusterSize Minimum cluster size.
#' @importFrom gplots heatmap.2
#' @return Heatmap of topic results.
#' @export
topic_heatmap <- function(df, method.dist = "euclidean", method.tree = "ward.D2", minClusterSize = 30){
  stopifnot(class(df) == "data.frame" | class(df) == "matrix" | class(df) == "topics")
  if (class(df) == "topics"){
    df <- df$omega
  }

  mycol = colorRampPalette(c("dark blue", "cyan", "yellow", "red"))(256)
  clusters = topic_clusters(df, method.dist = method.dist, method.tree = method.tree, minClusterSize = minClusterSize)
  #colors = c("#114477", "#4477AA", "#77AADD", "#117755", "#44AA88", "#99CCBB", "#777711", "#AAAA44", "#DDDD77", "#771111", "#AA4444", "#DD7777", "#771144", "#AA4477", "#DD77AA")
  colors = c("royalblue3",
             "mediumpurple4",
             "navajowhite2",
             "chocolate",
             "firebrick",
             "yellow2",
             "aquamarine",
             "orange1",
             "olivedrab2",
             "darkgreen",
             "pink",
             "black",
             "navy",
             "khaki3",
             "lightsteelblue1")
  clusters.col = colors[clusters]
  heatmap.2(
    t(df),
    key = TRUE,
    key.xlab = 'Identity',
    key.ylab = '',
    xlab = "feature",
    ylab = "topic",
    density.info='none',
    scale = 'none',
    trace = 'none',
    symbreaks = F,
    revC = F,
    cexRow = 0.7,
    cexCol = 0.35,
    symkey = 0,
    dendrogram = 'column',
    hclustfun = function(m)hclust(m, method = 'ward.D2'),
    distfun = function(m)dist(m,method = 'euclidean'),
    col = mycol,
    Rowv = F,
    ColSideColors = clusters.col)
}


#' Compute cluster matrix.
#'
#' This function is used to pool clustered features by adding the gene expression values within each cluster.
#' @param object Expression data.frame, matrix or object of class spaceST.
#' @param clusters Integer vector specifying cluster identity of each feature.
#' @return Integer vector specifying cluster identity of each feature.
#' @export
cluster_matrix <- function(object, clusters) {
  UseMethod("cluster_matrix")
}
#' @export
cluster_matrix.default <- function(object, clusters){
  if (!(class(object) %in% c("data.frame", "matrix"))){
    stop("Wrong input format.")
  }
  clust.matrix <- rowsum(t(object), clusters)
  clust.matrix <- as.data.frame(t(clust.matrix))
  colnames(clust.matrix) <- paste("c", 1:length(unique(clusters)), sep = "")
  return(clust.matrix)
}
#' @export
cluster_matrix.spaceST <- function(object, clusters){
  stopifnot(class(object) == "spaceST")
  clust.matrix <- cluster_matrix.default(object@expr, ifelse(is.null(clusters), object@meta.data$clusters, clusters))
}


#' Extract top features
#'
#' @param object Object of class spaceST.
ExtractTopFeatures <- function(object,
                               top_features = 1000,
                               method = "poisson",
                               options = "min",
                               shared = FALSE) {
  UseMethod("ExtractTopFeatures")
}
#' @export
ExtractTopFeatures.spaceST <- function(object,
                                       top_features = 1000,
                                       method = "poisson",
                                       options = "min",
                                       shared = FALSE) {
  top.features <- ExtractTopFeatures(object@lda.results$theta,
                                     top_features = top_features,
                                     method = method,
                                     options = options,
                                     shared = shared)
  return(top.features)
}
