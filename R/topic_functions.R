#' Helper function to compute lda using cellTree package.
#'
#' This helper function is used to compute topics for an expression dataset using lda modelling from cellTree package.
#' @name LDA
#' @param object Object of class spaceST.
#' @param min.topic,max.topic Integer values specifying min and max number of topics in the "maptpx" model.
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
#' @rdname LDA
#' @export
topic_compute.default <- function(
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
#' @rdname LDA
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
    d <- as.matrix(object@norm.data)
  } else if (datatype == "expr") {
    d <- as.matrix(object@expr)
  }
  lda.results <- topic_compute.default(
    d,
    min.topic = 2,
    max.topic = 15,
    num.genes = NULL,
    method = "maptpx",
    sd.filter = F,
    log.scale = F
  )
  object@lda.results <- lda.results
  object <- clusterST(object = object, minClusterSize = minClusterSize)
  object@meta.data$minClusterSize <- minClusterSize
  return(object)
}


#' Calculate clusters based on topics.
#'
#' This function is used to cluster features based on a topics matrix.
#' @param object Object of class spaceST.
#' @param method.dist Set distance method.
#' @param method.tree Set clustering method.
#' @param minClusterSize Integer value specifying the minimum cluster size allowed.
#' @param custom.matrix Provide a matrix to perform hierarchical clustering on.
#' @importFrom dynamicTreeCut cutreeDynamic
#' @return Integer vector specifying cluster identity of each feature.
#' @export
clusterST <- function(
  object = NULL,
  method.dist = "euclidean",
  method.tree = "ward.D2",
  minClusterSize = 30,
  custom.matrix = NULL
){
  stopifnot(!is.null(object) | !is.null(custom.matrix))
  if (!is.null(custom.matrix)) {
    omega <- custom.matrix
    my.dist = dist(omega, method = method.dist)
    my.tree = hclust(my.dist, method = method.tree)
    clusters = unname(cutreeDynamic(my.tree, distM = as.matrix(my.dist), verbose = 0, minClusterSize = minClusterSize))
    return(clusters)
  } else {
    omega <- object@lda.results$omega
    my.dist = dist(omega, method = method.dist)
    my.tree = hclust(my.dist, method = method.tree)
    clusters = unname(cutreeDynamic(my.tree, distM = as.matrix(my.dist), verbose = 0, minClusterSize = minClusterSize))
    object@meta.data$clusters <- clusters
    return(object)
  }
}


#' Plot heatmap of lda results.
#'
#' This function is used to plot a heatmap of lda results.
#' @param object Object of class spaceST.
#' @param method.dist Set distance method.
#' @param method.tree Set clustering method.
#' @param minClusterSize Minimum cluster size.
#' @param cols Set cluster colors.
#' @importFrom grDevices colorRampPalette
#' @importFrom stats hclust dist
#' @return Heatmap of topic results.
#' @export
topic_heatmap <- function(object, method.dist = "euclidean", method.tree = "ward.D2", minClusterSize = 30, cols = NULL){
  stopifnot(class(object) == "spaceST")
  df <- object@lda.results$omega
  mycol = colorRampPalette(c("dark blue", "cyan", "yellow", "red"))(256)
  clusters = clusterST(object, method.dist, method.tree, minClusterSize)@meta.data$clusters
  if (!is.null(cols)) {
    clusters.col <- cols
  } else {
    clusters.col <- brewer.pal(n = 9, name = "Set1")
  }
  clusters.col = clusters.col[clusters]

  gplots::heatmap.2(
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
#' @param top.features The top features in each cluster k that are selected based on the feature's ability
#' to distinguish cluster k from cluster 1, …, K for all cluster k < l. Default: 1000.
#' @param method The underlying model assumed for KL divergence measurement. Two choices considered are
#' "bernoulli" and "poisson". Default: poisson.
#' @param option if "min", for each cluster k, we select features that maximize the minimum KL divergence
#' of cluster k against all other clusters for each feature. If "max", we select features that maximize
#' the maximum KL divergence of cluster k against all other clusters for each feature.
#' @param shared if TRUE, then we report genes that can be highly expressed in more than one cluster.
#' Else, we stick to only those genes that are highest expressed only in a specific cluster.
#' @seealso \link[CountClust]{ExtractTopFeatures}
#' @export
ExtractTopFeaturesST <- function(
  object,
  top_features = 1000,
  method = "poisson",
  options = "min",
  shared = FALSE
 ) {
  theta <- object@lda.results$theta
  top.features <- CountClust::ExtractTopFeatures(theta,
                                     top_features = top_features,
                                     method = method,
                                     options = options,
                                     shared = shared)
  genes <- rownames(theta)
  res <- apply(top.features$indices, 1, function(x) {
    genes[x]
  })
  colnames(res) <- paste("Topic", 1:ncol(res), sep = "")
  return(res)
}
