#' spaceST network graphs
#'
#' @description Visualize SNN networks for spatial data
#' @param object Object of class spaceST.
#' @param subset_by Select attribute to subset igraph by, default NULL (options: "clusters", "replicate").
#' @param color_by Select attribute to color vertices by (options: "clusters", "replicate", "clust_and_rep).
#' @param mode Select mode for network object, default NULL (options: "spatial", "circle", "kamadakawai", ...
#' For more options, see available placement algorithms in the sna package).
#' @param select.rep.group Select replicate by replicate ID.
#' @param select.clust.group Select cluster by cluster number.
#' @param vertex.size Size of nodes/vertices.
#' @param edge.color Color edges of network.
#' @param edge.size Select size for edges in network.
#' @param clusters.snn Cluster data using network based approach.
#' @param cols character; Specify colors of grouping variables.
#' @param ... Parameters passed to ggnet2
#' @importFrom ggnet ggnet2
#' @importFrom igraph graph.adjacency E induced_subgraph V V<-
#' @importFrom graphics plot
#' @export
SpatialNetwork <- function(
  object,
  subset_by = NULL,
  color_by = "clusters",
  select.clust.group = 1,
  select.rep.group = NULL,
  mode = "fruchtermanreingold",
  vertex.size = 1,
  edge.color = "black",
  edge.size = 0.1,
  clusters.snn = F,
  cols = NULL,
  ...
) {
  if (is.null(object@snn)) {
    stop("No snn matrix present in spaceST object. Run CalcSNNspaceST() to compute SNN graph.")
  }
  stopifnot(
    !is.null(object@meta.data$clusters)
  )
  snn <- object@snn
  # Compute network
  net <- graph.adjacency(
    adjmatrix = as.matrix(snn),
    mode = "undirected",
    weighted = TRUE,
    diag = FALSE
  )
  if (clusters.snn) {
    V(net)$clusters <- object@meta.data$clusters.snn
  } else {
    V(net)$clusters <- object@meta.data$clusters
  }
  V(net)$replicates <- as.numeric(object@coordinates$replicate)
  V(net)$x <- object@coordinates$x
  V(net)$y <- 35 - object@coordinates$y
  if (mode == "tsne") {
    stopifnot(length(object@tsne) > 0)
    V(net)$x <- object@tsne[, 1]
    V(net)$y <- object@tsne[, 2]
  }
  # Subset graph
  if (!is.null(subset_by)) {
    if (subset_by == "clusters") {
      net <- induced_subgraph(net, V(net)[clusters %in% select.clust.group])
    } else if (subset_by == "replicate") {
      if (!select.rep.group %in% unique(object@coordinates$replicate)) {
        stop(paste("Not a valid replicate ID. Available IDs:", print(unique(object@coordinates$replicate))))
      }
      net <- induced_subgraph(net, V(net)[replicates == select.rep.group])
    } else if (subset_by == "clust_and_rep") {
      if (!select.rep.group %in% unique(object@coordinates$replicate)) {
        stop(paste("Not a valid replicate ID. Available IDs:", print(unique(object@coordinates$replicate))))
      }
      net <- induced_subgraph(net, V(net)[replicates == select.rep.group])
      net <- induced_subgraph(net, V(net)[clusters %in% select.clust.group])
    }
  }
  # Define layout
  if (!is.null(mode)) {
    if (mode == "spatial" | mode == "tsne") {
      net.layout <- c("x", "y")
    } else {
      net.layout <- mode
    }
  }
  # Define coloring
  if (color_by == "clusters") {
    col <- V(net)$clusters
  } else if (color_by == "replicate") {
    col <- V(net)$replicates
  } else {
    col <- color_by
  }
  if (is.null(cols)) {
    cols <- c("royalblue3", "mediumpurple4", "navajowhite2", "chocolate", "firebrick", "yellow2"," aquamarine", "orange1", "olivedrab2", "darkgreen", "pink", "black", "navy", "khaki3", "lightsteelblue1")
  }
  col <- cols[col]
  pnet <- ggnet2(asNetwork(net),
         mode = net.layout,
         color = col,
         size = vertex.size,
         edge.color = edge.color,
         edge.size = edge.size,
         #palette = cols,
         ...
         )
  plot(pnet)
}
