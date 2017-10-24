#' Interpolate object values over "cell population"
#'
#' @description Function used to interpolate values, e.g. gene expression values, xCell scores or topic proportions
#' over dot pattern.
#' @param val.df Data.frame with object values.
#' @param section.input Data.frame with dot pattern coordinates.
#' @param x Number of grids along x axis used for raster.
#' @param y Number of grids along y axis used for raster.
#' @importFrom akima interp
#' @return data.frame with input values for feature coordinates interpolated over 2D scatter.
#' @export
interpolate_2D_data <- function(val.df, section.input, x, y, iterval = NULL){
  x1 = as.numeric(val.df[,1])
  x1 = x1[!is.na(x1)]
  y1 = as.numeric(val.df[,2])
  y1 = y1[!is.na(y1)]
  w1 = as.numeric(val.df[,3])
  s1 =  interp(x1, y1, w1, nx = x, ny = y)
  mat.1 = s1$z
  mat.1 = as.data.frame(mat.1)
  colnames(mat.1) = c(1:y)
  mat.1 = mat.1[,rev(colnames(mat.1)),]
  mat.1 = as.matrix(mat.1)
  col = as.numeric(mat.1)
  if (!is.null(iterval)) {
    set = col[section.input[,4]]
    section.xyz.value = cbind(section.input[,1:2], rep(iterval, nrow(section.input)), set)
    df <- as.data.frame(section.xyz.value)
    colnames(df) <- c("x", "y", "z", "val")
  } else {
    set = col[section.input[,3]]
    section.xyz.value = cbind(section.input[,1:2], set)
    df <- as.data.frame(section.xyz.value)
    colnames(df) <- c("x", "y", "val")
  }
  return(df)
}
