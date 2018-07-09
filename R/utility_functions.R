#' Quality Control statistics
#'
#' @description Given a data.frame with gene expression values (genes as rows and observations as columns), this function will return a
#' data.frame with unique genes and transcript counts for each observation, as well as the sample id
#' @export
#' @param df Input gene expression data.frame
#' @return Data.frame with unique genes per feature and transcripts per feature
ST_statistics <- function(df){
  unique.genes.per.feature <- apply(df, 2, function(x) sum(x > 0))
  transcripts.per.feature <- colSums(df)
  return(data.frame(unique.genes.per.feature, transcripts.per.feature))
}

#' Quality Control plot
#'
#' @description This function plots histograms of unique genes per feature and transcripts per featrure.
#' @export
#' @param df Input quality control data.frame.
#' @param separate Should replicates be seperated?
#' @param bins Select number of bins.
#' @param ... Parameters passed to geom_histogram.
#' @importFrom ggplot2 ggplot aes_string geom_histogram facet_grid
#' @return Plot of histograms.
ST_statistics_plot <- function(df, separate = F, bins = NULL, ...){
  if (!is.null(bins)) {
    ranges <- apply(df[, 1:2], 2, range)
    binwidths <- round(c(ranges[2, 1] - ranges[1, 1], ranges[2, 2] - ranges[1, 2])/bins)
  }
  p.list <- lapply(c("unique.genes.per.feature", "transcripts.per.feature"), function(type) {
    if (!is.null(bins)) {
      bw <- ifelse(type == "unique.genes.per.feature", binwidths[1], binwidths[2])
    } else {
      bw <- NULL
    }
    p <- ggplot(df, aes_string(type)) +
      geom_histogram(color = "black", fill = ifelse(type == "unique.genes.per.feature", scales::muted("red"), scales::muted("blue")), binwidth = bw, ...)
    if (separate) {
      p <- p + facet_grid(~samples)
    }
    return(p)
  })
  plot(cowplot::plot_grid(plotlist = p.list, nrow = 2))
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
get_coordinates <- function(x, delimiter = "_") {
  UseMethod("get_coordinates")
}

#' @rdname coordinates
#' @export
get_coordinates.default <- function(x, delimiter = "_"){
  coords <- do.call(rbind, strsplit(colnames(x), split = delimiter))
  colnames(coords) <- c("replicate", "x", "y")
  coords <- data.frame(coords, stringsAsFactors = F)
  coords$x <- as.numeric(coords$x)
  coords$y <- as.numeric(coords$y)
  return(coords)
}

#' @rdname coordinates
#' @export
get_coordinates.spaceST <- function(x, delimiter = "_"){
  x@coordinates
}

#' @rdname coordinates
#' @export
get_coordinates.dgCMatrix <- function(x, delimiter = "_"){
  get_coordinates.default(x, delimiter = delimiter)
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
#' @param delimiter Delimiter used in header.
#' @return Merged and filtered dataframe.
merge_exp_list <- function(x,
                  unique.genes = 0,
                  min.exp = 0,
                  min.features = 0,
                  filter.genes = NULL,
                  delimiter = "_"){
  if (class(x) == "list") {
    stopifnot(length(x) > 1)
    df.A <- as.matrix(x[[1]])

    # Check headers
    coords <- do.call(rbind, strsplit(colnames(df.A), split = delimiter))
    if (dim(coords)[2] == 2) {
      xy <- paste(1, colnames(df.A), sep = delimiter)
    } else if (dim(coords)[2] == 3) {
      xy <- colnames(df.A)
    }

    for (i in 2:length(x)){
      df.B <- as.matrix(x[[i]])
      coords <- do.call(rbind, strsplit(colnames(df.B), split = delimiter))
      if (dim(coords)[2] == 2) {
        xy <- c(xy, paste(i, colnames(df.B), sep = delimiter))
      } else if (dim(coords)[2] == 3) {
        xy <- c(xy, colnames(df.B))
      }
      df.A <- data.frame(merge(df.A, df.B, by = "row.names", all = TRUE), row.names = 1)
    }
    df.A[is.na(df.A)] <- 0
    colnames(df.A) <- xy
    all.samples.matrix <- as(as.matrix(df.A), "dgCMatrix")
  } else {
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

  # Filter out low quality spots
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
#' Duplicated gene names will be aggregated.
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
ensembl2hgnc.list <- function(df, organism = "human") {
  exp.list <- list()
  for (i in 1:length(df)) {
    exp.list[[i]] <- ensembl2hgnc.default(df[[i]], organism = organism)
  }
  return(exp.list)
}

#' @rdname convert
#' @export
ensembl2hgnc.spaceST <- function(spST) {
  spST@expr <- as(as.matrix(ensembl2hgnc.default(spST@expr)), "dgCMatrix")
  return(spST)
}

#' Calculate cp10k
#'
calc_cp10k <- function(df) {
  norm.factors <- colSums(df)
  norm.data <- t(t(df)/norm.factors*1e4)
}

