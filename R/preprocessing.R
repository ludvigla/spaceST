#' Initiialize and setup the spaceST object
#'
#' Initialized the spaceST object from an expression matrix with genes as rows and array coordinates as columns. When the object is initialized,
#' some optional filtering settings can be applied to remove genes and/or spots of low quality.
#
#' @param raw.data A list of expression matrices/data.frames with genes as rows and feature coordinates as columns.
#' Headers should be formatted with a unique symbol, x, and y coordinate seperated by a delimiter.
#' @param unique.genes The lowest number of unique genes allowed in a feature.
#' @param min.exp Integer value specifying the lowest expression value allowed at least min.features number of features.
#' @param min.features Integer value specifying the lowest number of features with at least min.exp.
#' @param filter.genes A character vector specifying genes that should be filtered from the expression data.
#' @param delimiter Delimiter specifying header format.
#' @return A spaceST object containing a merged expression matrix of gene expression counts.
#' @importFrom CountClust BatchCorrectedCounts
#' @export
CreatespaceSTobject <- function(
  raw.data,
  unique.genes = 0,
  min.exp = 0,
  min.features = 0,
  filter.genes = NULL,
  delimiter = "_"
  ) {
  spaceST.version <- packageVersion("spaceST")
  if (!(class(raw.data) %in% c("list", "data.frame", "matrix"))) {
    return("Wrong input format")
  }
  object <- new(
    Class = "spaceST",
    filter.settings = list(
      unique.genes = unique.genes,
      min.exp = min.exp,
      min.features = min.features,
      filter.genes = filter.genes),
    version = spaceST.version
  )
  object@expr <- as(merge_exp_list(
    raw.data, unique.genes,
    min.exp, min.features,
    filter.genes = filter.genes
  ), "dgCMatrix")
  object@coordinates <- get_coordinates(object@expr, delimiter)
  object@status.expr <- "filtered raw data"
  return(object)
}
