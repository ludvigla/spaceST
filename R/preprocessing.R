#' Filter ST expression data.
#'
#' @description This function is used to filter ST expression data used for 3D analysis. Features with low number of unique
#' genes and genes with low expression values will be removed. By default, ribosomal proteins and MALAT1 genes are also removed
#' @param raw.data A list of expression matrices/data.frames with genes as rows and feature coordinates as columns
#' Required structure of column names is: "replicate number_x_y"
#' @param unique.genes The lowest number of unique genes allowed in a feature
#' @param min.exp Integer value specifying the lowest expression value allowed at least min.features number of features (default 2)
#' @param min.features Integer value specifying the lowest number of features with at least min.exp (default 15)
#' @param filter.genes A character vector specifying genes that should be filtered from the expression data
#' @param batch.correct Logical specifying if data should be batch corrected
#' @return A spaceST object containing filtered gene expression data for multiple replicates
#' @importFrom CountClust BatchCorrectedCounts
#' @export
CreatespaceSTobject <- function(
  raw.data,
  unique.genes = 500,
  min.exp = 2,
  min.features = 15,
  filter.genes = NULL
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
  object@coordinates <- get_coordinates(object@expr)
  object@status.expr <- "filtered raw data"
  return(object)
}
