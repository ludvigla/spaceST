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
#' @param bg.black Set background color to black.
#' @param xlim,ylim Set limits of x/y axes. [default: xlim = c(1, 34), ylim = c(1, 36)]
#' @param arrange Arrange plots.
#' @param col.title Give the color legend a title.
#' @param scale Set this parameter to "colwise" to scale each column separately.
#' @param size Set size of markers.
#' @param hide.dropouts Logical; Set to TRUE if you want to hide spots with 0 values.
#' @param ncols Number of columns in arranged plot table.
#' @param ... Parameters passed to geom_point.
#' @importFrom grid rasterGrob
#' @importFrom jpeg readJPEG
#' @importFrom ggplot2 ggplot geom_point annotation_custom theme element_blank labs theme_void scale_color_gradientn element_text
#' @importFrom grid unit
#' @importFrom grDevices rgb
#' @importFrom stats na.omit
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
  bg.black = FALSE,
  xlim = c(1, 33),
  ylim = c(1, 35),
  ncols = NULL,
  arrange = TRUE,
  col.title = NULL,
  scale = "",
  size = 3,
  hide.dropouts = F,
  ...
) {
  stopifnot(class(object) == "spaceST")
  if (length(value) == ncol(object@expr)) {
    val = value
  } else {
    if (type == "expr") {
      val = object@expr[value, ]
      if (hide.dropouts) {
        val[val == 0] <- NA
        print(sum(is.na(val)))
      }
    } else if (type == "norm.data") {
      val = object@norm.data[value, ]
      if (hide.dropouts) {
        val[val == 0] <- NA
        print(sum(is.na(val)))
      }
    } else if (type == "pca") {
      val = object@reducedDims$x[, value]
    } else if (type == "topics") {
      val = object@lda.results$omega[, value]
    } else if (dim(value)[2] %in% c(2, 3)) {
      d <- dim(value)[2]
      if (scale == "colwise") {
        val <- apply(value, 2, function(x) {
         return((x - min(x))/(max(x) - min(x)))
        })
      } else {
        val <- value
        x <- as.vector(val)
        val <- matrix((x - min(x))/(max(x) - min(x)), ncol = d)
      }
      if (d == 2) {
        val <- cbind(val, 0)
      }
      val <- rgb(val, maxColorValue = 1)
    } else {
      stop("type ", type, " not valid. Select one of 'expr', 'norm.data' or 'pca'")
    }
  }
  # Combine array coordinates with value
  gg.df <- data.frame(object@coordinates, val)
  max.col <- max(val, na.rm = T)
  min.col <- min(val, na.rm = T)
  pal <- palette.select(palette)
  reps <- unique(gg.df$replicate)

  # Convert jpegs to grobs
  if (!is.null(HE.list)) {
    if (length(HE.list) != length(reps)) {
      stop("Number of images in HE.list does not match number of samples.")
    } else {
      grobs.list <- lapply(1:length(HE.list), function(i) {
        HE_img <- jpeg::readJPEG(HE.list[[i]])
        g <- rasterGrob(HE_img, width = unit(1, "npc"), height = unit(1, "npc"), interpolate = TRUE)
        return(g)
      })
    }
  }

  # Plot spatial heatmap
  p.list <- lapply(1:length(reps), function(i) {
    subset.df <- na.omit(subset(gg.df, replicate == reps[i]))
    cols <- rgb(pal(seq(0, 1, length.out = 10)), maxColorValue = 255)
    if (invert.heatmap) {
      cols <- rev(cols)
    }
    if (class(val) == "character") {
      p <- ggplot(subset.df, aes(x, 36 - y))
    } else {
      p <- ggplot(subset.df, aes(x, 36 - y, color = val))
    }
    p <- p +
      scale_x_continuous(limits = xlim, expand = c(0, 0)) +
      scale_y_continuous(limits = ylim, expand = c(0, 0))
    if (!is.null(HE.list)) {
      p <- p + annotation_custom(grobs.list[[i]], -Inf, Inf, -Inf, Inf)
    }
    if (nrow(subset.df) > 1007) {
      size = 0.1
    }
    if (class(val) == "character") {
      p <- p + geom_point(color = subset.df$val, size = size, ...) +
        theme_void()
    } else {
      p <- p + geom_point(size = size, ...) +
        theme_void() +
        labs(color = ifelse(!is.null(col.title), col.title, value))
        if (scale == "colwise") {
          p <- p + scale_color_gradientn(colours = cols)
        } else {
          p <- p +scale_color_gradientn(colours = cols, limits = c(min.col, max.col))
        }
    }
    if (hide.legend) {
      p <- p + guides(color = FALSE)
    }
    if (bg.black) {
      p <- p + theme(plot.background = element_rect(fill = "black"),
                     legend.text = element_text(colour = "white"),
                     legend.title = element_text(colour = "white"))
    }
    return(p)
  })
  if (arrange) {
    cowplot::plot_grid(plotlist = p.list, ncol = ncols)
  } else {
    return(p.list)
  }
}


#' Visualize spatial distribution of clusters.
#'
#' @param object Object of class spaceST.
#' @param clusters Integer/numeric vector specifying clusters.
#' @param HE.list List of paths to HE images in jpeg format that should be used as a background for the
#' spatial heatmap.
#' @param arrange Arrange plots.
#' @param ncols Number of columns in arranged plot table.
#' @param xlim,ylim Set limits of x/y axes. [default: xlim = c(1, 34), ylim = c(1, 36)]
#' @param cols Set custom cluster colors by specifying a character vector with color codes. You can specify
#' specific colors for each cluster by naming the character vector with the cluster ids.
#' @param ... additional parameters passed to geom_point(), defining the attributes of the feature coordinate points.
#' @importFrom grid rasterGrob
#' @importFrom ggplot2 ggplot geom_point annotation_custom theme element_blank labs scale_color_gradientn theme_void scale_color_manual
#' @return Plot of clustered featured overlayed on black and white HE image.
#' @export
spatial.clusters <- function(object,
                             clusters = NULL,
                             HE.list = NULL,
                             arrange = T,
                             ncols = NULL,
                             xlim = c(1, 33),
                             ylim = c(1, 35),
                             cols = NULL,
                             ...) {
  # Combine array coordinates with value
  reps <- unique(object@coordinates$replicate)
  if (!is.null(clusters)) {
    gg.df <- data.frame(object@coordinates, cluster = clusters)
  } else if (!length(object@meta.data$clusters) == 0) {
    gg.df <- data.frame(object@coordinates, cluster = object@meta.data$clusters)
  } else {
    stop("No cluster vector present in spaceST object. Please run topic_compute or provide clusters ...")
  }

  # Convert jpegs to grobs
  if (!is.null(HE.list)) {
    if (length(HE.list) != length(reps)) {
      stop("Number of images in HE.list does not match number of samples.")
    } else {
      grobs.list <- lapply(1:length(HE.list), function(i) {
        HE_img <- jpeg::readJPEG(HE.list[[i]])
        g <- rasterGrob(HE_img, width = unit(1, "npc"), height = unit(1, "npc"), interpolate = TRUE)
        return(g)
      })
    }
  }

  # Plot spatial heatmap
  p.list <- lapply(1:length(reps), function(i) {
    subset.df <- na.omit(subset(gg.df, replicate == reps[i]))
    p <- ggplot(subset.df, aes(x, 36 - y, color = factor(cluster)))
    if (!is.null(HE.list)) {
      p <- p + annotation_custom(grobs.list[[i]], -Inf, Inf, -Inf, Inf)
    }
    p <- p + geom_point(...) +
      theme_void() +
      labs(color = "cluster") +
      scale_x_continuous(limits = xlim, expand = c(0, 0)) +
      scale_y_continuous(limits = ylim, expand = c(0, 0))
    if (!is.null(cols) & length(unique(gg.df$cluster)) <= length(cols)) {
      p <- p + scale_color_manual(values = cols)
    }

    return(p)
  })
  if (arrange) {
    cowplot::plot_grid(plotlist = p.list, ncol = ncols)
  } else {
    return(p.list)
  }
}
