#' Generate binary dot pattern
#'
#' @description Given a black and white HE image, this function can be used to create a scatter representatio of cell nuclei ditributed
#' over the tissue image. The image needs to be completely free from noise surrounding the tissue area and make sure that the colorspace
#' is set to grey.
#' @importFrom jpeg readJPEG
#' @importFrom ggplot2 ggplot aes theme geom_point element_blank
#' @param img File path or list of file paths for BW image(s).
#' @param limit Defines cutoff on the greyscale for what to interpret as a valid coordinate.
#' @param rownum Number of rows in output grid, i.e. the maximum number of points to keep.
#' @examples
#' library(spaceST)
#' # Create scatter pattern and plot results
#'
#' scatter <- scatter_HE(img = <path to black and white HE image>, show.plot = T)
#' @return Data.frame with x, y coordinates representing cell nuclei on an HE image.
#' @export
scatter_HE <- function(img, coords = NULL, limit = 0.5, rownum = 5e4, show.plot = FALSE, offset_x = 0.75, offset_y = 0.75) {
  UseMethod("scatter_HE")
}
#' @export
scatter_HE.default <- function(img,
                               coords = NULL,
                               limit = 0.5,
                               rownum = 5e4,
                               show.plot = FALSE,
                               offset_x = 0.75,
                               offset_y = 0.75){
  if (class(img) != "character"){
    stop("Wrong input format.")
  } else if (!file.exists(img)) {
    stop("No such file exists.")
  }
  bw.image = jpeg::readJPEG(img)
  img = which(bw.image < limit , arr.ind = TRUE)

  test.y = attributes(bw.image)$dim[1] / 34
  test.x = attributes(bw.image)$dim[2] / 32

  V1 = (img[,2] / test.x) + 1
  V2 = (img[,1] / test.y) + 1

  #if (!is.null(coords)) {
  #  V1 = scale2range(V1, (min(coords[,1]) - offset_x), (max(coords[, 1]) + offset_x))
  #  V2 = scale2range(V2, (min(coords[,2]) - offset_y), (max(coords[, 2]) + offset_y))
  #}

  img = cbind(V1, V2)
  set.seed(1)
  if (rownum < nrow(img)) {
    img <- img[sample(1:nrow(img), size = rownum, replace = FALSE), ]
  }
  img <- as.data.frame(img)
  colnames(img) <- c("x", "y")
  set.seed(NULL, normal.kind = "default")
  if (show.plot) {
    p <- ggplot(img, aes(x, -y)) +
      geom_point(stat = "identity", size = 0.5) +
      theme(axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.x = element_blank(),
            panel.background = element_blank(),
            panel.grid = element_blank(),
            panel.grid.major = element_blank())
    plot(p)
  }
  return(img)
}
#' @export
scatter_HE.list <- function(img, limit = 0.5, rownum = 5e4, show.plot = FALSE, offset){
  grid_list <- list()
  for (i in 1:length(img)) {
    grid_list[[i]] <- scatter_HE(img[[i]], limit = limit, rownum = rownum, offset = offset)
  }
  return(grid_list)
}

#' rasterize scatter data from HE image
#'
#' @description This function is used to rasterize scatter data from a black and white HE image generated with function `scatter_HE()`.
#' @param scatter Data.frame with binary dot pattern.
#' @param coords Coordinate data.frame with x, y values. Use `get_coordinates()` function to obtain feature coordinates for a
#' specific replicate.
#' @importFrom raster cellFromXY raster
#' @importFrom sp coordinates
#' @return data.frame with x, y coordinates and grid cell numbers. Column names for x, y coordinates gives the dimensions of the grid.
#' @export
rasterize_scatter <- function(scatter,
                              coords,
                              plot = FALSE, raster_offset = 0.5) {
  UseMethod("rasterize_scatter")
}
#' @export
rasterize_scatter.default <- function(scatter, coords){
  stopifnot(is.data.frame(scatter) && nrow(scatter) != 0,
            is.data.frame(coords) && nrow(coords) != 0)
  max.x = max(coords[,1]) + 0.5
  min.x = min(coords[,1]) - 0.5
  max.y = max(coords[,2]) + 0.5
  min.y = min(coords[,2]) - 0.5
  x = round((max(coords[,1]) - min(coords[,1]) + 1), digits = 0)*4
  res = (max.x - min.x) / x
  r = raster(xmn = min.x, ymn = min.y, xmx = max.x, ymx = max.y, res = res)
  r[] = 0

  tab = table(cellFromXY(r, scatter))

  r[as.numeric(names(tab))] = tab

  pixel.centers = coordinates(r)
  set1 = scatter[,1:2]
  set2 = pixel.centers[,1:2]

  new.set2 <- apply(set1, 1, function(x) which.min(colSums((t(set2) - x)^2)))

  section.input = cbind(scatter[,1], scatter[,2], new.set2)
  y = r@nrows
  colnames(section.input) = c(x, y, 'grid.cell')
  return(as.data.frame(section.input))
}
#' @export
rasterize_scatter.list <- function(scatter, coords) {
  stopifnot(is.list(scatter),
            is.list(coords))
  raster_list <- list()
  for (i in 1:length(scatter)) {
    raster_list[[i]] <- rasterize_scatter(scatter[[i]], coords[[i]])
  }
  return(raster_list)
}
