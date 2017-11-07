#' Filter ST expression data.
#'
#' @description This function is used to filter ST expression data used for 3D analysis. Features with low number of unique
#' genes and genes with low expression values will be removed. By default, ribosomal proteins and MALAT1 genes are also removed
#' @param raw.data A list of expression matrices/data.frames with genes as rows and feature coordinates as columns
#' Required structure of column names is: "replicate number_x_y"
#' @param unique.genes The lowest number of unique genes allowed in a feature
#' @param min.exp Integer value specifying the lowest expression value allowed at least min.features number of features (default 2)
#' @param min.features Integer value specifying the lowest number of features with at least min.exp (default 15)
#' @param filter.data A character vector specifying genes that should be filtered from the expression data (default RP)
#' @param batch.correct Logical specifying if data should be batch corrected
#' @return A spaceST object containing filtered gene expression data for multiple replicates
#' @importFrom CountClust BatchCorrectedCounts
#' @examples
#' library(STanalysis3D)
#' data(bcST)
#'
#' # Filter with default dettings
#' ST.object <- get.spaceST(bcST)
#' ST.object
#'
#' # Filter data to keep features with more than 300 unique genes and genes
#'  with at least 1 count in 10 features
#' # and disable filtering of ribosomal protein coding genes/MALAT1.
#' ST.object <- get.spaceST(bcST,
#'                          unique.genes = 300,
#'                          min.exp = 1,
#'                          min.features = 10,
#'                          filter.data = NULL)
#' ST.object
#' @export
CreatespaceSTobject <- function(
  raw.data,
  unique.genes = 500,
  min.exp = 2,
  min.features = 15,
  filter.data = RP,
  batch.correct = FALSE
  ) {
  spaceST.version <- packageVersion("spaceST")
  if (!(class(raw.data) %in% c("list", "data.frame", "matrix"))) {
    return("Wrong input format")
  }
  object <- new(
    Class = "spaceST",
    raw.data = raw.data,
    filter.settings = list(
      unique.genes = unique.genes,
      min.exp = min.exp,
      min.features = min.features,
      filter.data = filter.data),
    version = spaceST.version
  )
  object.raw.data <- object@raw.data
  object@expr <- merge_exp_list(
    raw.data,unique.genes,
    min.exp, min.features,
    filter.data = filter.data
  )
  if (batch.correct){
    object@corrected <- as.data.frame(
      t(BatchCorrectedCounts(t(object@expr),
                             get_coordinates(object@expr)[, 1],
                             use_parallel = T)))
  }
  object@coordinates <- get_coordinates(object@expr)
  return(object)
}
