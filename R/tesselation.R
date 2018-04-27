#' Voronoi tesselation plot
#'
#' @description Compute voronoi tesselation and visualize ST data as a colored mosaic
#' @param r data.frame with grouped x, y values used with geom_polygon
#' @param datapolyM convex hull polygon used to mask area outside tissue
#' @param outline.color set color of polygon outline
#' @param bg background color
#' @param xlim set limits of x axis
#' @param ylim set limits of y axis
#' @importFrom ggplot2 ggplot geom_polygon scale_x_continuous scale_y_continuous aes theme_void theme element_rect margin
tesselation.plot <- function(
  r,
  datapolyM = NULL,
  outline.color = NA,
  bg = "black",
  xlim = c(0, 34),
  ylim = c(0, 36),
  alpha = 1
) {
  if (is.na(outline.color)) {
    outline.color <- r$value
  }
  p <- ggplot(r, aes(x = x, y = y)) +
    geom_polygon(aes(group = id), fill = r$value, color = outline.color, alpha = alpha) +
    scale_x_continuous(limits = c(xlim[1], xlim[2]), expand = c(0, 0)) +
    scale_y_continuous(limits = c(ylim[1], ylim[2]), expand = c(0, 0)) +
    theme_void() +
    theme(plot.margin = margin(0,0,-0.15,-0.15, "cm"),
          plot.background = element_rect(fill = bg)) +
    guides(fill = FALSE)
  if (!is.null(datapolyM)) {
    p <- p + geom_polygon(data = datapolyM, aes(x = x, y = ylim[2] - y), fill = bg, alpha = 1)
  }
  return(p)
}


#' Convert tiles and add color.
#'
#' @description convert tile.list object produced with the deldir package into a data.frame structure that is
#' interpretable by geom_polygon.
#' @param tiles tile.list object produced with the tile.list function from teh deldir package.
#' @param v vector of colors representing some value to be visualized.
#' @param xlim,ylim set limits of x/y axes.
tiles2polygon <- function(tiles, v, xlim, ylim) {
  r <- do.call(rbind, lapply(1:length(tiles), function(i) data.frame(x = tiles[[i]]$x, y = tiles[[i]]$y, value = rep(v[i], length(tiles[[i]]$x)), id = rep(i, length(tiles[[i]]$x)), stringsAsFactors = F)))
  pts <- as.data.frame(do.call(rbind, lapply(1:length(tiles), function(i) tiles[[i]]$pt)))
  r$x[r$x < xlim[1]] <- xlim[1]
  r$x[r$x > xlim[2]] <- xlim[2]
  r$y[r$y < ylim[1]] <- ylim[1]
  r$y[r$y > ylim[2]] <- ylim[2]
  return(r)
}

#' Plot tesselated object.
#'
#' @description Sisualize some value on an ST grid using voronoi tesselation.
#' @param object Sbject of class spaceST.
#' @param select.columns Select columns of target matrix [default = 1:3].
#' @param xlim Set limits of x axis [default = c(0, 36)].
#' @param ylim Set limits of y axis [default = c(0, 36)].
#' @param scale.by.column Should columns be scaled independently? [default = FALSE].
#' @param datatype Select target matrix from spaceST object [default = "topics"].
#' @param add.points Should points be added to the plot?
#' @param outline.color Set polygon outline color [default = NA].
#' @param bg Select background color [default = "black"].
#' @param inv.col Invert colors [default = FALSE].
#' @param return.list Return list of plots instead of combining them into a grid [default = FALSE].
#' @param cell.col Color of points (spots or scatter) [default = "cyan"].
#' @param cell.size Size of points (spots or scatter) [default = 0.01].
#' @param cell.alpha Color of points (spots or scatter) [default = 0.1].
#' @param fill.empty Fill empty array spots with 0 values.
#' @importFrom deldir deldir tile.list
#' @importFrom grDevices chull
#' @importFrom ggplot2 geom_point aes
#' @importFrom magrittr %>%
#' @export
tessViz <- function(
  object,
  select.columns = 1:3,
  xlim = c(0, 34),
  ylim = c(0, 36),
  scale.by.column = FALSE,
  datatype = "topics",
  add.points = F,
  outline.color = NA,
  bg = "black",
  inv.col = F,
  alpha = 1,
  scatter = NULL,
  return.list = FALSE,
  cell.col = "cyan",
  cell.size = 0.01,
  cell.alpha = 0.1,
  palette = "offwhite.to.black",
  max.num.cols = 20,
  polygon.mask = TRUE,
  fill.empty = FALSE
) {
  stopifnot(class(object) == "spaceST",
            datatype %in% c("topics", "reducedDims", "t-SNE"))
  all.coords <- object@coordinates
  if (datatype == "topics") {
    stopifnot(length(object@lda.results$omega) > 0)
    all.z = object@lda.results$omega
  } else if (datatype == "reducedDims") {
    stopifnot(length(object@reducedDims) > 0)
    all.z = object@reducedDims
  } else if (datatype == "t-SNE") {
    stopifnot(length(object@tsne) > 0)
    all.z <- object@tsne
  }

  # Generate empty marix to fill empty spot positions
  xy <- apply(all.coords[, 2:3], 2, round)
  all_xy <- c()
  for (i in 1:33) {
    for (j in 1:35) {
      all_xy <- c(all_xy, c(i, j))
    }
  }
  all_xy <- t(matrix(all_xy, nrow = 2))

  # Empty list to store plots in
  plot.list <- list()
  for (n in unique(all.coords$replicate)) {
    # Subset coordinates per replicate
    ind <- which(all.coords$replicate == n)
    coords <- all.coords[ind, 2:3]
    x <- coords[, 1]
    y <- coords[, 2]
    y <- 36 - y

    # subset z values per replicate
    z <- all.z[ind, ]

    # Scale data
    if (scale.by.column) {
      z <- apply(z, 2, scale2range, a = 0, b = 1)
    } else {
      z <- matrix(scale2range(as.numeric(z), a = 0, b = 1), ncol = ncol(z))
    }

    # Create polygon from convex hull
    if (polygon.mask) {
      tissue.polygon <- coords[chull(x, y),] %>% apply(2, rev)
      tissue.polygon <- rbind(tissue.polygon, tissue.polygon[1, ])
      rim  <- c(xlim[1],ylim[1], xlim[2],ylim[1], xlim[2],ylim[2], xlim[1],ylim[2], xlim[1],ylim[1]) %>% matrix(ncol = 2, byrow = T)
      datapolyM <- rbind(rim, tissue.polygon) %>% as.data.frame()
      names(datapolyM) <- c("x","y")
      datapolyM$id <- "a"
    } else {
      datapolyM <- NULL
    }

    # Add 0 values to empty spots
    if (fill.empty) {
      all_xy_subset <- all_xy[!apply(all_xy, 1, paste, collapse = "x") %in% paste(round(x), round(y), sep = "x"), ]
      x <- c(x, all_xy_subset[, 1])
      y <- c(y, all_xy_subset[, 2])
      z <- rbind(z, matrix(0, nrow = nrow(all_xy_subset), ncol = ncol(z)))
    }

    # Run voronoi tesselation
    vtess <- suppressWarnings(suppressMessages(deldir(x, y)))
    tiles <- tile.list(vtess)

    # Plot each column in Grayscale if select columns is not specified
    if (!is.null(select.columns)) {
      z <- z[, select.columns]
      if (length(select.columns) > 3) {
        warning("More than 3 columns selected. Using first 3 columns ...")
      }
      if (ncol(z) < 3) {
        while (ncol(z) < 3) {
          z <- cbind(z, 0)
        }
      }
      if (inv.col) {
        z <- apply(z, 2, function(x) 1 - x)
      }

      # Convert
      v <- rgb(z)
      r <- tiles2polygon(tiles = tiles, v = v, xlim = xlim, ylim = ylim)
      p <- tesselation.plot(r = r,
                            datapolyM = datapolyM,
                            outline.color = outline.color,
                            bg = bg, xlim = xlim, ylim = ylim, alpha = alpha)

      if (add.points) {
        if (!is.null(scatter)) {
          p <- p + geom_point(data = scatter[[n]], aes(x, y), color = cell.col, size = cell.size, alpha = cell.alpha)
        } else {
          p <- p + geom_point(data = coords, aes(x, 36 - y), size = cell.size, color = cell.col, alpha = cell.alpha)
        }
      }
      inside_plot <- p
    } else {
      inside_plot <- list()
      for (i in 1:ncol(z)) {
        if (i > max.num.cols) {
          warning("Maxumim number of plots exceeded ...")
          break
        }
        v <- z[, i]
        if (inv.col) {
          v <- 1 - v
        }

        # Convert
        pal <- palette.select(palette)
        v <- rgb(pal(v), maxColorValue = 255)
        #v <- grey(v)
        r <- tiles2polygon(tiles = tiles, v = v, xlim = xlim, ylim = ylim)
        p <- tesselation.plot(r = r,
                              datapolyM = datapolyM,
                              outline.color = outline.color,
                              bg = bg, xlim = xlim, ylim = ylim, alpha = alpha)

        if (add.points) {
          if (!is.null(scatter)) {
            p <- p + geom_point(data = scatter[[n]], aes(x, y), color = cell.col, size = cell.size, alpha = cell.alpha)
          } else {
            p <- p + geom_point(data = coords, aes(x, 36 - y), size = cell.size, alpha = cell.alpha)
          }
        }
        inside_plot[[i]] <- p
      }
    }
    if ("list" %in% class(inside_plot)) {
      inside_plot <- cowplot::plot_grid(plotlist = inside_plot)
    }
  plot.list[[n]] <- inside_plot
  }
  if (length(plot.list) == 1) {
    plot.list[[1]]
  } else if(return.list) {
    return(plot.list)
  } else {
    cowplot::plot_grid(plotlist = plot.list)
  }
}


# Scale values
scale2range <- function(x, a, b){
  res <- (b-a)*(x - min(x))/(max(x) - min(x)) + a
  return(res)
}

#' Create palettes
#'
#' @description returns a palette using the coloRamp function
#' @param palette select palette [options: "GnBu", "the.cols", "Spectral", "offwhite.to.black"]
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRamp
#' @importFrom viridis inferno magma plasma
#' @export
palette.select <- function(palette) {
  palettes <- list(
    green.to.blue = colorRamp(brewer.pal(9,"GnBu")),
    the.cols = colorRamp(c(rgb(255,255,217, maxColorValue=255),
                           rgb(65,182,196, maxColorValue=255),
                           rgb(8, 29, 88, maxColorValue=255)),
                         space="Lab"),
    spectral = colorRamp(brewer.pal(9,"Spectral")),
    offwhite.to.black = colorRamp(c(rgb(220,220,220, maxColorValue=255),
                                    rgb(0, 0, 0, maxColorValue=255)),
                                  space="Lab"),
    viridis = colorRamp(viridis(9)),
    cividis = colorRamp(cividis(9)),
    magma = colorRamp(magma(9)),
    plasma = colorRamp(plasma(9))
  )
  return(palettes[[palette]])
}
