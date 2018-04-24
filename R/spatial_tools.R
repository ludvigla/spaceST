#' Visualize spatial heatmaps
#'
#' Generates a scatter plot of array spots colored by some value.
#' @param object Object of class spaceST.
#' @param value Target vector to visualize. This value can be chosen from any slot of the spaceST object
#' containing data linked to array spots.
#' @param type Select dataset where the value can be found [options: "expr", "norm.data", "pca"]
#' @param HE.list List of paths to HE images in jpeg format that should be used as a background for the
#' spatial heatmap.
#' @param palette Color palette used for spatial heamtap [options: ""expr", "green.to.blue", "spectral", "offwhite.to.black", "
#' "viridis", "magma", "plasma", "cividis"]
#' @param invert.heatmap Invert color gradient.
#' @param hide.legend Exclude legend.
#' @param ... Parameters passed to geom_point.
#' @importFrom grid rasterGrob
#' @importFrom jpeg readJPEG
#' @importFrom ggplot2 ggplot geom_point annotation_custom theme element_blank
#' @rdname heatmaps
#' @export
spatial.heatmap <- function(
  object,
  value,
  type = "expr",
  HE.list = NULL,
  palette = "spectral",
  invert.heatmap = FALSE,
  hide.legend = FALSE,
  ...
) {
  stopifnot(class(object) == "spaceST")
  if (type == "expr") {
    val = object@expr[value, ]
  } else if (type == "norm.data") {
    val = object@norm.data[value, ]
  } else if (type == "pca") {
    val = object@reducedDims[, value]
  } else {
    stop("type ", type, " not valid. Select one of 'expr', 'norm.data' or 'pca'")
  }

  # Combine array coordinates with value
  gg.df <- data.frame(object@coordinates, val)
  pal <- palette.select(palette)
  reps <- unique(gg.df$replicate)

  # Convert jpegs to grobs
  if (!is.null(HE.list)) {
    if (length(HE.list) != length(reps)) {
      stop("Number of images in HE.list does not match number of samples.")
    } else {
      grobs.list <- lapply(HE.list, function(x) {
        img <- jpeg::readJPEG(HE.list[[i]])
        g <- rasterGrob(HE_img, width = unit(1, "npc"), height = unit(1, "npc"), interpolate = TRUE)
        return(g)
      })
    }
  }

  # Plot spatial heatmap
  p.list <- lapply(1:length(reps), function(i) {
    subset.df <- subset(gg.df, replicate == reps[i])
    cols <- rgb(pal(seq(0, 1, length.out = 10)), maxColorValue = 255)
    if (invert.heatmap) {
      cols <- rev(cols)
    }
    p <- ggplot(subset.df, aes(x, 36 - y, color = val))
    if (!is.null(HE.list)) {
      p <- p + annotation_custom(grobs.list[[i]], -Inf, Inf, -Inf, Inf)
    }
    p <- p + geom_point(...) +
      theme_void() +
      labs(color = value) +
      scale_color_gradientn(colours = cols)
    if (hide.legend) {
      p <- p + guides(color = FALSE)
    }
    return(p)
  })
  cowplot::plot_grid(plotlist = p.list)
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
#' @export
spatial.clusters <- function(clusters, coords, HE, size = 4, type = "ST_grid", ...) {
  stopifnot(is.integer(clusters) | is.numeric(clusters),
            is.data.frame(coords),
            is.character(HE),
            file.exists(HE))
  clust.coords <- cbind(clusters, coords)
  if (type == "ST_grid") {
    p <- ggplot(clust.coords, aes(x = x-1, y = 35-y, color = factor(clusters, levels = (1:max(unique(clusters))))))
  } else if (type == "pixels") {
    p <- ggplot(clust.coords, aes(x = x, y = y, color = factor(clusters, levels = (1:max(unique(clusters))))))
  }
  if (!is.null(HE)) {
    HE_img = jpeg::readJPEG(HE)
    g <- rasterGrob(HE_img, width = unit(1, "npc"), height = unit(1, "npc"), interpolate = TRUE)
    p <- p + annotation_custom(g, -Inf, Inf, -Inf, Inf) +
      geom_point(size = size, ...)
    if (type == "pixels") {
      print(dim(HE_img)[2]*2.5)
      p <- p + scale_x_continuous(limits = c(0, dim(HE_img)[2]), expand = c(0, 0)) +
        scale_y_continuous(limits = c(0, dim(HE_img)[1]), expand = c(0, 0))
    }
  } else {
    p <- p + geom_point(size = size, ...)
  }
  p <- p + theme(axis.line = element_blank(),
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
    labs(color = "unique genes (mean)")
  if (type == "ST_grid") {
    p <- p + scale_color_manual(values = c("1" = "royalblue3",
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
                                  "15" = "lightsteelblue1"))
  } else {
    p <- p + scale_color_brewer(palette = "Set1")
  }
  plot(p)
}
