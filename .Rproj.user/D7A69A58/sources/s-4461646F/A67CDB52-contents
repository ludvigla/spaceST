#' Visualize spatial heatmaps
#'
#' @description Visualize spatial features on top of HE image.
#' @param x Vector with numerical values with the same length as coords. Values could represent gene expresion values,
#' topic proportions, xCell scores etc.
#' @param coords data.frame with feature x, y coordinates.
#' @param raster data.frame with x, y coordinates and grid cell numbers produced with the scatter_rasterize() function.
#' @param HE path to HE image
#' @param alpha set alpha level of heatmap
#' @param size set size of points
#' @param mirror_x flip along x axis
#' @param mirror_y flip along y axis
#' @importFrom grid rasterGrob
#' @importFrom ggplot2 ggplot geom_point annotation_custom theme element_blank
#' @return Heatmap of spatially distributed scores.
#' @export
spatial.heatmap <- function(x,
                            coords,
                            raster,
                            HE = NULL,
                            alpha = 0.5,
                            size = 1,
                            mirror_x = F,
                            mirror_y = F,
                            title = NULL,
                            shift_x = -1,
                            shift_y = 35,
                            lim_x = c(0, 32),
                            lim_y = c(0, 34),
                            max = 1) {
  stopifnot(is.numeric(x) | is.integer(x),
            is.data.frame(coords),
            is.data.frame(raster))
  if (is.character(HE)) {
    stopifnot(file.exists(HE))
  }
  score.df <- cbind(coords, x)
  colnames(score.df) <- c("x", "y", "val")
  score.df <- interpolate_2D_data(score.df, raster, x = as.numeric(colnames(raster)[1]), y = as.numeric(colnames(raster)[2]))
  p <- ggplot(score.df, aes(x = x+shift_x, y = shift_y-y, color = val))
  if (!is.null(HE)) {
    HE_img = jpeg::readJPEG(HE)
    g <- rasterGrob(HE_img, width = unit(1, "npc"), height = unit(1, "npc"), interpolate = TRUE)
    p <- p + annotation_custom(g, -Inf, Inf, -Inf, Inf) +
      geom_point(alpha = alpha, size = size)
  } else {
    p <- p + geom_point(alpha = alpha, size = size)
  }
  if (class(score.df$val) == "integer") {
    p <- p + scale_x_continuous(limits = lim_x, expand = c(0, 0)) +
      scale_y_continuous(limits = lim_y, expand = c(0, 0)) +
      scale_color_brewer(palette = "set1") +
      ggplot2::theme(axis.line = ggplot2::element_blank(),
                     panel.grid.major = ggplot2::element_blank(),
                     panel.grid.minor = ggplot2::element_blank(),
                     panel.border = ggplot2::element_blank(),
                     panel.background = ggplot2::element_blank(),
                     axis.title.x = element_blank(),
                     axis.text.x = element_blank(),
                     axis.ticks.x = element_blank(),
                     axis.title.y = element_blank(),
                     axis.text.y = element_blank(),
                     axis.ticks.y = element_blank())
    plot(p)
  } else if (class(score.df$val) == "numeric") {
    p <- p + scale_x_continuous(limits = lim_x, expand = c(0, 0)) +
      scale_y_continuous(limits = lim_y, expand = c(0, 0)) +
      ggplot2::theme(axis.line = ggplot2::element_blank(),
                     panel.grid.major = ggplot2::element_blank(),
                     panel.grid.minor = ggplot2::element_blank(),
                     panel.border = ggplot2::element_blank(),
                     panel.background = ggplot2::element_blank(),
                     axis.title.x = element_blank(),
                     axis.text.x = element_blank(),
                     axis.ticks.x = element_blank(),
                     axis.title.y = element_blank(),
                     axis.text.y = element_blank(),
                     axis.ticks.y = element_blank()) +
      scale_colour_gradientn(colours = c("red","#E69F00","#F0E68C","cyan","darkblue"),
                             values=c(1, 0.8, 0.6, 0.4,0.2,0),
                             limits = c(0, max),
                             guide = "colorbar",
                             na.value = "white") +
      ggtitle(title)
    if (mirror_x) {
      p <- p + scale_x_reverse()
    }
    if (mirror_y) {
      p <- p + scale_y_reverse()
    }
    plot(p)
  }
}

#' Visualize spatial distribution of clusters
#'
#' @description Given a black and white HE image, aligned feature coordinates and a cluster vector, this function
#' can be used to visualize the spatial distribution of clusters.
#'
#' NOTE: The input image should be cropped to the centers of the frame corners, i.e. the outermost layer of features.
#' This is a requirement for the
#' @param clusters Integer vector specifying feature cluster identity for corresponding feature coordinates.
#' @param coords data.frame with x, y values for feature coordinates.
#' @param HE Path to cropped black and white HE image.
#' @param size numerical specifying coordinate size.
#' @param ... additional parameters passed to geom_point(), defining the attributes of the feature coordinate points.
#' @importFrom grid rasterGrob
#' @importFrom ggplot2 ggplot geom_point annotation_custom theme element_blank
#' @return Plot of clustered featured overlayed on black and white HE image.
#' @examples
#' library(spaceST)
#'
#' data(bcST)
#' spaceST.object <- CreatespaceSTobject(bcST, corrected = T)
#'
#' # Run factor analysis and cluster features
#' spaceST.object <- topic_compute(spaceST.object)
#' clusters <- topic_clusters(topic.df)
#'
#' # Get coordinates from spaceST object
#' coords <- get_coordinates(spaceST)[, 1:2]
#'
#' # Plot spatial heatmap
#' spatial.clusters(clusters, coords, HE = "/path/to/HE")
#' @export
spatial.clusters <- function(clusters, coords, HE, size = 4, ...) {
  stopifnot(is.integer(clusters),
            is.data.frame(coords),
            is.character(HE),
            file.exists(HE))
  clust.coords <- cbind(clusters, coords)
  p <- ggplot(clust.coords, aes(x = x-1, y = 35-y, color = factor(clusters, levels = (1:max(unique(clusters))))))
  if (!is.null(HE)) {
    HE_img = jpeg::readJPEG(HE)
    g <- rasterGrob(HE_img, width = unit(1, "npc"), height = unit(1, "npc"), interpolate = TRUE)
    p <- p + annotation_custom(g, -Inf, Inf, -Inf, Inf) +
      geom_point(size = size, ...)
  } else {
    p <- p + geom_point(size = size, ...)
  }
  p <- p + scale_x_continuous(limits = c(0, 32), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 34), expand = c(0, 0)) +
    scale_color_manual(values = c("1" = "royalblue3",
                                  "2" = "mediumpurple4",
                                  "3" = "navajowhite2",
                                  "4" = "chocolate",
                                  "5" = "firebrick",
                                  "6" = "yellow2",
                                  "7" = "aquamarine",
                                  "8" = "orange1",
                                  "9" = "olivedrab2",
                                  "10" = "darkgreen",
                                  "11" = "pink",
                                  "12" = "black",
                                  "13" = "navy",
                                  "14" = "khaki3",
                                  "15" = "lightsteelblue1")) +
    theme(axis.line = element_blank(),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.border = element_blank(),
                   panel.background = element_blank(),
                   axis.title.x = element_blank(),
                   axis.text.x = element_blank(),
                   axis.ticks.x = element_blank(),
                   axis.title.y = element_blank(),
                   axis.text.y = element_blank(),
                   axis.ticks.y = element_blank()) +
    labs(color = "clusters")
  plot(p)
}
