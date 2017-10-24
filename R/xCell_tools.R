#' Calculate xCell scores for cluster matrix
#'
#' @description This function is used to subset xCell results based on celltype.
#' @param clust.df xCell data.frame with clusters or features as columns and celltypes as rows.
#' @return A subset data.frame with the celltypes of interest.
#' @export
subset.xCell <- function(clust.df, type = "lymphoids") {
    celltype <- get(type)
    return(clust.df[which(rownames(clust.df) %in% celltype), ])
}

#' Rearrange xCell output from clusters to features
#'
#' @description This function can rearrange the clusters in the xCell score matrix back into features.
#' During the rearrangement, the cluster score is assigned to each feature within that cluster. This assumes
#' that every feature within each cluster has the same celltype mixture, which is a rough approximation.
#' However, high xCell scores indicate a strong presence of a certain celltype in that cluster even though
#' the exact celltype positioning is unknown.
#' @param clusters Integer vector specifying cluster positions.
#' @param features Character vector with feature names of the original data.frame that was used for pooling.
#' @return A rearranged xCell matrix with features as columns and genes as rows.
#' @export
clust2feature <- function(xCell.matrix, clusters, features) {
  if (!(class(xCell.matrix) %in% c("data.frame", "matrix"))) {
    stop("Wrong input format of xCell.matrix")
  }
  if (!(class(clusters) %in% c("numeric", "integer"))) {
    stop("Wrong input format of clusters")
  }
  if (class(features) != "character") {
    stop("Wrong input format of features")
  }
  celltypes <- rownames(xCell.matrix)
  # Obtain sorting vector
  sorted <- c()
  for (c in sort(unique(clusters))){
    positions <- (1:length(clusters))[clusters == c]
    names(positions) <- rep(as.character(c), length(positions))
    sorted <- c(sorted, positions)
  }

  # Now, we need to divide the xCell score of each cluster and recreate the original matrix design
  # with features as columns. First, we divide the clusters and bind in an ordered sequence.
  clust.matrix.list <- list()

  n <- 1
  for (c in sort(unique(clusters))){
    clustnum <- sum(clusters == c)
    divided.matrix <- matrix(rep(xCell.matrix[, paste("c", c, sep = "")], clustnum), ncol = clustnum) # add /clustnum?
    colnames(divided.matrix) <- rep(c, clustnum)
    clust.matrix.list[[n]] <- divided.matrix
    n <- n + 1
  }

  ordered.matrix <- do.call(cbind, clust.matrix.list)
  arranged.matrix <- matrix(nrow = nrow(ordered.matrix), ncol = ncol(ordered.matrix))

  # Rearrange matrix to match original format
  for (i in 1:length(sorted)){
    arranged.matrix[, sorted[i]] <- ordered.matrix[, i]
  }

  rownames(arranged.matrix) <- celltypes
  colnames(arranged.matrix) <- features

  return(arranged.matrix)
}

#' xCell heatmap
#'
#' @description Given an xCell output data.frame, this function can be used to generate a heatmap of xCell scores for celltypes
#' and clusters. You can subset the data using the default groups provided in the STanalysis3D package or provide a character vector
#' of cell type names. Clusters and cell types are automatically ordered with dendrograms using hclust (method = "Ward.D2").
#' @param df data.frame with xCell output, i.e. cell types as rows and clusters as columns.
#' @param group Used to subset xCell data. Preset groups are; lymphoids, stem, stromal, myeloids, scores, others.
#' @param scale If TRUE, xCell scores will be scaled to a [0, 1] scale.
#' @param compare If TRUE, xCell scores will be scaled among all celltypes regardless of subtype group, making results
#' comparable between subtype groups.
#' @param col2 Set "low" and "high" colors for heatmap.
#' @importFrom ggdendro dendro_data segment
#' @importFrom ggplot2 element_blank theme ggplot geom_tile aes coord_flip scale_x_discrete scale_y_discrete
#' unit scale_fill_gradient2 element_rect labs ggplotGrob
#' @importFrom gtable gtable_add_rows gtable_add_cols gtable_add_grob
#' @importFrom grid grid.newpage grid.draw
#' @importFrom reshape2 melt
#' @return Heatmap of xCell scores.
#' @export
xCell.heatmap <- function(df, group = NULL, scale = FALSE, compare = F, col.names = NULL, col2 = c("#4B0082", "#F0E442")) {
  if (!is.null(col.names)) {
    colnames(df) <- col.names
  }
  if (!is.null(group)) {
    df <- df[rownames(df) %in% group, ]
  }
  x <- as.matrix(df)
  dd.col <- as.dendrogram(hclust(dist(t(x))))
  col.ord <- order.dendrogram(dd.col)
  dd.row <- as.dendrogram(hclust(dist(x)))
  row.ord <- order.dendrogram(dd.row)
  xx <- x[row.ord, col.ord]
  xx_names_1 <- colnames(xx)
  xx_names_2 <- rownames(xx)
  df.dendro <- as.data.frame(xx)
  colnames(df.dendro) <- xx_names_1
  df.dendro$cells <- xx_names_2
  df.dendro$cells <- with(df.dendro, factor(cells, levels = cells, ordered = TRUE))
  mdf <- melt(df.dendro, id.vars = "cells")
  if (scale) {
    mdf$value <- scale2range(mdf$value, a = 0, b = 1)
  }
  ddata_x <- dendro_data(dd.row)
  ddata_y <- dendro_data(dd.col)

  theme_none <- theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.title.x = element_text(colour=NA),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.line = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank())

  p1 <- ggplot(mdf, aes(variable, cells)) +
    geom_tile(aes(fill = value), colour = "white") +
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(x = "cluster", y = "celltype", fill = "score") +
    theme(legend.position = "left")
  if (compare) {
    p1 <- p1 + scale_fill_gradient2(low = col2[1], high = col2[2], midpoint = (max(x)/2))
  } else {
    p1 <- p1 + scale_fill_gradient2(low = col2[1], high = col2[2], midpoint = (max(mdf$value)/2))
  }
  g <- ggplotGrob(p1)
  p2 <- ggplot(segment(ddata_x)) +
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
    theme_none + coord_flip() +
    theme(plot.margin = unit(c(0, 0, 0, -0.3), "cm"))
  g2 <- ggplotGrob(p2)
  g <- gtable_add_cols(x = g, unit(1.5, "in"))
  g <- gtable_add_grob(x = g, grobs = g2, t = 1, l = 7, b = 8, r = -1)
  p3 <- ggplot(segment(ddata_y)) +
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
    theme_none +
    theme(plot.margin = unit(c(0, 0.3, -0.7, 0), "cm"))
  g3 <- ggplotGrob(p3)
  g <- gtable_add_rows(g, unit(1, "in"), 0)
  g <- gtable_add_grob(x = g, grobs = g3, t = 1, b = 6, l = 6, r = 6)
  grid.newpage()
  grid.draw(g)
}

