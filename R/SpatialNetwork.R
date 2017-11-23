#' spaceST network graphs
#'
#' @description Visualize SNN networks for spatial data
#' @param subset_by Select attribute to subset igraph by, default NULL (options: "clusters", "replicate")
#' @param color_by Select attribute to color vertices by (options: "clusters", "replicate", "clust_and_rep)
#' @param mode Select mode for network object, default NULL (options: "spatial", "circle", "kamadakawai", ...
#' For more options, see available placement algorithms in the sna package)
#' @param select.rep.group Select replicate by replicate ID
#' @param select.clust.group Select cluster by cluster number
#' @param vertex.size Size of nodes/vertices
#' @param edge.color Color edges of network
#' @param edge.size Select size for edges in network
#' @param ... Parameters passed to ggnet2
#' @importFrom ggnet ggnet2
#' @importFrom igraph graph.adjacency V E induced_subgraph
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
  cols <- c("royalblue3", "mediumpurple4", "navajowhite2", "chocolate", "firebrick", "yellow2","aquamarine", "orange1", "olivedrab2", "darkgreen", "pink", "navy", "khaki3", "lightsteelblue1")
  if (length(unique(col)) > 14) {
    cols = c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788")
  }
  col <- cols[col]
  pnet <- ggnet2(net,
         mode = net.layout,
         color = col,
         size = vertex.size,
         edge.color = edge.color,
         edge.size = edge.size,
         #palette = cols,
         ...
         )
  if (mode == "spatial") {
    pnet <- pnet + scale_x_continuous(limits = c(0, 32)) +
      scale_y_continuous(limits = c(0, 35))
  }
  plot(pnet)
}
