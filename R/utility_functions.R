#' Quality Control statistics
#'
#' @description Given a data.frame with gene expression values (genes as rows and observations as columns), this function will return a
#' data.frame with unique genes and transcript counts for each observation, as well as the sample id
#' @export
#' @param df Input gene expression data.frame
#' @return Data.frame with unique genes per feature and transcripts per feature
ST_statistics<- function(df){
  stopifnot(
    is.data.frame(df)
  )
  samples <- do.call(rbind, strsplit(colnames(df), split = "_"))[, 1]
  unique.genes.per.feature <- apply(df, 2, function(x) sum(x > 0))
  transcripts.per.feature <- colSums(df)
  return(data.frame(unique.genes.per.feature, transcripts.per.feature, samples))
}

#' Quality Control plot
#'
#' @description This function plots histograms of unique genes per feature and transcripts per featrure.
#' @export
#' @param df Input quality control data.frame.
#' @param separate Should replicates be treated separately?
#' @return Plot of histograms.
#' @importFrom cowplot plot_grid
#' @importFrom ggplot2 element_blank theme aes geom_histogram facet_grid scale_color_manual geom_vline element_line
#' @importFrom scales muted
ST_statistics.plot <- function(df, separate = F){
  theme_none <- theme(axis.line = element_line(colour = "black"),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.border = element_blank(),
                      panel.background = element_blank())
  if (separate == T) {
    p1 <- ggplot(df, aes(x = unique.genes.per.feature)) +
      geom_histogram(bins = 30, fill = muted("dark blue"), color = "black") +
      theme_none +
      scale_color_manual(name = "", values = c(median = "black")) +
      facet_grid(~samples)

    p2 <- ggplot(df, aes(x = transcripts.per.feature)) +
      geom_histogram(bins = 30, fill = muted("red"), color = "black") +
      theme_none +
      scale_color_manual(name = "", values = c(median = "black")) +
      facet_grid(~samples)
  } else {
    p1 <- ggplot(df, ggplot2::aes(x = unique.genes.per.feature)) +
      geom_histogram(bins = 30, fill = muted("dark blue"), color = "black") +
      geom_vline(aes(xintercept = median(df$unique.genes.per.feature), colour = "median"),
                          linetype="dashed",
                          size=1) +
      theme_none +
      ggplot2::scale_color_manual(name = "", values = c(median = "black"))

    p2 <- ggplot(df, aes(x = transcripts.per.feature)) +
      geom_histogram(bins = 30, fill = muted("red"), color = "black") +
      geom_vline(aes(xintercept = median(df$transcripts.per.feature), colour = "median"),
                          linetype="dashed",
                          size=1) +
      theme_none +
      scale_color_manual(name = "", values = c(median = "black"))
  }
  p <- plot_grid(p1, p2, nrow = 2)
  plot(p)
}

#' Obtain coordinates from expression data list of expression data.frame
#'
#' @title get_coordinates: obtain spatial coordinates
#' @description This function is used to collect all coordinates from a list of expression datasets or an expression dataset
#' @export
#' @rdname coordinates
#' @param x List of expression data.frames, data.frame or matrix.
#' @param rep Get replicate column
#' @return Data frame with replicate numbers, x and y coordinates
get_coordinates <- function(x, rep = TRUE) {
  UseMethod("get_coordinates")
}

#' @rdname coordinates
#' @export
get_coordinates.matrix <- function(x, rep = TRUE){
  get_coordinates.data.frame(x, rep = TRUE)
}

#' @rdname coordinates
#' @export
get_coordinates.data.frame <- function(x, rep = TRUE){
    ids <- colnames(x)
    coords <- as.data.frame(do.call(rbind, strsplit(ids, split = "_")))
    coords[, 2] <- as.numeric(as.character(coords[, 2]))
    coords[, 3] <- as.numeric(as.character(coords[, 3]))
    colnames(coords) <- c("replicate", "x", "y")
    if (rep) {
      return(coords)
    } else {
      return(coords[, 2:3])
    }
}

#' @rdname coordinates
#' @export
get_coordinates.spaceST <- function(x, rep = TRUE){
  if (nrow(x@corrected) == 0) {
    get_coordinates.data.frame(x@expr, rep = rep)
  } else {
    get_coordinates.data.frame(x@corrected, rep = rep)
  }
}

#' @rdname coordinates
#' @export
get_coordinates.list <- function(x, rep = TRUE) {
  coords_list <- list()
  for (i in 1:length(x)) {
    coords_list[[i]] <- get_coordinates(x[[i]], rep = rep)
  }
  return(as.data.frame(do.call(rbind, coords_list)))
}

#' Merge and filter expression datasets
#'
#' @description This function is used to merge datasets in a list of expression datasets.
#' Genes which are not present in some datasets will get 0 values in those datasets. The merged dataset is
#' subsequently filtered from features with few unique genes and genes with low expression counts.
#' @export
#' @rdname merge
#' @param x List of exrpression data frames, data.frame or matrix.
#' @param unique.genes Integer value specifying number of unique genes allowed in a feature.
#' @param min.exp Integer value specifying lowest expression value allowed in min.features number of features.
#' @param min.features Integer value specifying number of features allowed with min.exp count.
#' @param filter.genes Character vector specifying genes that should be filtered out.
#' @return Merged and filtered dataframe.
merge_exp_list <- function(x,
                  unique.genes = 0,
                  min.exp = 0,
                  min.features = 0,
                  filter.genes = NULL){
  if (class(x) == "list") {
    df.A <- x[[1]]
    for (i in 2:length(x)){
      df.B <- x[[i]]
      df.A <- data.frame(merge(df.A, df.B, by = "row.names", all = TRUE), row.names = 1)
    }
    df.A[is.na(df.A)] <- 0
    all.samples.matrix <- as(as.matrix(df.A), "dgCMatrix")
  } else if (class(x) %in% c("data.frame", "matrix")) {
    all.samples.matrix <- as(as.matrix(x), "dgCMatrix")
  }

  if (!is.null(filter.genes)){
    removed.genes <- -grep(rownames(all.samples.matrix), pattern = filter.genes, perl = TRUE)
    if(length(removed.genes) > 0) {
      all.samples.matrix = all.samples.matrix[-grep(rownames(all.samples.matrix), pattern = filter.genes, perl = TRUE),]
    }
  }

  # Filter out low quality genes
  all.samples.matrix = all.samples.matrix[Matrix::rowSums(all.samples.matrix >= min.exp) >= min.features, ]

  # Filter out low quality features
  indices <- which(apply(all.samples.matrix, 2, function(x) sum(x > 0)) < unique.genes)
  if (length(indices) > 0){
    all.samples.matrix <- all.samples.matrix[, -indices]
  }
  return(all.samples.matrix)
}

#' Collect data from bioMmRt
#'
#' @description function used to collect data from biomaRt and save to data.frame
#' @export
#' @param gene.list character vector of gene ids
#' @param organism select organism (default hsapiens)
#' @param filters select filters from input data.frame
#' @param attributes select attributes to collect from biomaRt
collect_gene_data <- function(gene.list, organism = "hsapiens", filters = "ensembl_gene_id_version", attributes = c("ensembl_gene_id_version", "gene_biotype")) {
  mart <- character()
  class(mart) <- "try-error"
  trials <- 0
  while (class(mart) == "try-error") {
    mart <- try(useDataset(paste(organism, "_gene_ensembl", sep = ""), useMart("ensembl", ensemblRedirect = FALSE)), silent = TRUE)
    trials <- trials + 1
    if (trials == 5) {
      stop("Failed to connect to BiomaRt")
    }
  }
  G_list <- getBM(filters = filters, attributes = attributes, values = gene.list, mart = mart)
  if (nrow(G_list) == 0) {
    stop("No matches found. Check if the corect organism was selected.")
    }
  return(G_list)
}

#' Convert from ENSEMBL ID to HGNC symbol or MGI symbol
#'
#' @description Function used to convert ENSEMBL gene ids of an expression matrix into HGNC/MGI symbols.
#' @export
#' @rdname convert
#' @param df Data.frame or matrix with ENSEMBL ids as rownames.
#' @param organism Select organism database (human, mouse)
#' @return Data.frame or matrix with converted names.
#' @importFrom biomaRt getBM useDataset useMart
ensembl2hgnc <- function(df, organism = "human") {
  UseMethod("ensembl2hgnc")
}

#' @rdname convert
#' @export
ensembl2hgnc.default <- function(df, organism = "human") {
  mart <- character()
  class(mart) <- "try-error"
  trials <- 0
  if (organism == "human") {
    while (class(mart) == "try-error") {
      mart <- try(useDataset("hsapiens_gene_ensembl", useMart("ensembl")), silent = TRUE)
      trials <- trials + 1
      if (trials == 15) {
        stop("Failed to connect to BiomaRt")
      }
    }
    G_list <- getBM(filters= "ensembl_gene_id_version", attributes = c("ensembl_gene_id_version", "hgnc_symbol"), values = rownames(df), mart = mart)
    if (nrow(G_list) == 0) {
      stop("No matches found. Check if the corect organism was selected.")
    }
    G_list$hgnc_symbol[G_list$hgnc_symbol == ""] <- G_list$ensembl_gene_id_version[G_list$hgnc_symbol == ""]
  } else if(organism == "mouse") {
    while (class(mart) == "try-error") {
      mart <- try(useDataset("mmusculus_gene_ensembl", useMart("ensembl")), silent = TRUE)
      trials <- trials + 1
      if (trials == 15) {
        stop("Failed to connect to BiomaRt")
      }
    }
    G_list <- getBM(filters= "ensembl_gene_id_version", attributes = c("ensembl_gene_id_version", "mgi_symbol"), values = rownames(df), mart = mart)
    if (nrow(G_list) == 0) {
      stop("No matches found. Check if the corect organism was selected.")
    }
    G_list$mgi_symbol[G_list$mgi_symbol == ""] <- G_list$ensembl_gene_id_version[G_list$mgi_symbol == ""]
  }

  G_vec <- G_list[, 2]
  names(G_vec) <- G_list[, 1]
  df <- as.matrix(df)
  rownames(df)[!is.na(G_vec[rownames(df)])] <- G_vec[rownames(df)[!is.na(G_vec[rownames(df)])]]
  if (sum(duplicated(rownames(df))) > 0) {
    df <- rowsum(df, rownames(df))
  }
  return(as.data.frame(df))
}

#' @rdname convert
#' @export
ensembl2hgnc.list <- function(df) {
  exp.list <- list()
  for (i in 1:length(df)) {
    exp.list[[i]] <- ensembl2hgnc.default(df[[i]])
  }
  return(exp.list)
}

#' @rdname convert
#' @export
ensembl2hgnc.spaceST <- function(spST) {
  spST@expr <- ensembl2hgnc.default(spST@expr)
  return(spST)
}

#' Cast merged data to list of matrixes
#'
#' @description Cast any merged data.frame/matrix into a list of matrices for each replicate.
#' @export
#' @param df Merged data.frame/matrix with spaceST headers, i.e. features as columns.
#' @return List of matrices.
cast2list <- function(df) {
  samples <- as.integer(as.character(get_coordinates(df)[, 1]))
  exp.list <- list()
  for (i in unique(samples)) {
    indices <- (1:ncol(df))[which(samples == i)]
    exp.list[[i]] <- df[, indices]
  }
  return(exp.list)
}


#' Normalize using CPTK method
#'
#' @description Function used to normalize ST data using Counts Per Ten thousand method
#' @export
#' @rdname norm
#' @param expr_mat Expression data.frame or matrix.
#' @return Normalized data.frame or matrix.
calc_cpm <- function (expr_mat) {
  norm_factor <- colSums(expr_mat)
  return(t(t(expr_mat)/norm_factor)*10^4)
}
