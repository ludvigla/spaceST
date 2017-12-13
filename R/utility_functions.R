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
#' @param filter.data Character vector specifying genes that should be filtered out.
#' @return Merged and filtered dataframe.
merge_exp_list <- function(x,
                       unique.genes,
                       min.exp,
                       min.features,
                       filter.data){
  if (class(x) == "list") {
    for (i in 1:length(x)){
      df <- x[[i]]
      # Order gene names alphabetically
      df <- df[order(rownames(df)), ]
      # Add index row
      df <- as.data.frame(cbind(id = rownames(df), df))
      x[[i]] <- df
    }
    df.A <- x[[1]]
    for (i in 2:length(x)){
      df.B <- x[[i]]
      df.A <- merge(df.A, df.B, by = "id", all = TRUE)
    }
    gene.names <- df.A$id
    # remove id column
    sample.df <- df.A[, -1]
    sample.df <- apply(sample.df, 2, as.numeric)
    rownames(sample.df) <- gene.names
    sample.df[is.na(sample.df)] <- 0
  } else if (class(x) %in% c("data.frame", "matrix")) {
    sample.df <- as.data.frame(x)
  }
  if (!is.null(filter.data)){
    sample.df = sample.df[!rownames(sample.df) %in% filter.data,]
  }

  # Filter out low quality genes
  sample.df = sample.df[rowSums(sample.df >= min.exp) >= min.features,]
  sample.df <- as.data.frame(sample.df)

  # Filter out low quality features
  indices <- which(apply(sample.df, 2, function(x) sum(x > 0)) < unique.genes)
  if (length(indices) > 0){
    filtered.sample.df <- sample.df[, -indices]
  } else {
    filtered.sample.df <- sample.df
  }
  return(filtered.sample.df)
}

#'  pca of expression data
#'
#' @description This function is used to  the first two principal components of an expression dataset, or a comparison
#' between raw expression and corrected expression data.
#' @param df1 Expression dataset 1
#' @param df2 Expression dataset 2
#' @param samples Vector of replicate numbers of length ncol(filtered)
#' @param ... arguments passed to plotPCA function from scater
#' @importFrom scater newSCESet plotPCA
#' @importFrom cowplot plot_grid
#' @importFrom Biobase pData
#' @return Scatter plot of PC1 and PC2.
pca_plot <- function(df1, df2 = NULL, samples, ...){
  sce.raw = newSCESet(countData = df1)
  pData(sce.raw)$replicate = samples
  p1 <- plotPCA(sce.raw, pca_data_input = "counts", colour_by = "replicate", ...)
  if (!is.null(df2)) {
    sce.batch.corr = newSCESet(countData = df2)
    pData(sce.batch.corr)$replicate = samples
    p2 <- plotPCA(sce.batch.corr, pca_data_input="counts", colour_by = "replicate", ...)
    p <- plot_grid(p1, p2, nrow = 2, label_size = 10,
                            labels = c(paste("      raw data"),
                                       paste("batch corrected")))
  } else {
    p <- p1
  }
  plot(p)
}

#' Convert from ENSEMBL ID to HGNC symbol
#'
#' @description  Function used to convert ENSEMBL gene ids of an expression matrix into HGNC symbols.
#' @export
#' @rdname convert
#' @param df Data.frame or matrix with ENSEMBL ids as rownames.
#' @return Data.frame or matrix with converted names.
#' @importFrom biomaRt useDataset getBM useMart
ensembl2hgnc <- function(df) {
  UseMethod("ensembl2hgnc")
}

#' @rdname convert
#' @export
ensembl2hgnc.data.frame <- function(df) {
  mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id", "hgnc_symbol"), values = rownames(df), mart = mart)
  df <- cbind(ensembl_gene_id = rownames(df), df)
  merged.df <- merge(G_list, df, by = "ensembl_gene_id", all = T)
  merged.df$hgnc_symbol[which(merged.df$hgnc_symbol == "" | is.na(merged.df$hgnc_symbol))] <- merged.df$ensembl_gene_id[which(merged.df$hgnc_symbol == "" | is.na(merged.df$hgnc_symbol))]
  merged.df <- data.frame(aggregate(x = apply(merged.df[, 3:ncol(merged.df)], 2, function(x) as.numeric(as.character(x))), by = list(merged.df$hgnc_symbol), sum, na.rm = TRUE), row.names = 1)
  return(merged.df)
}

#' @rdname convert
#' @export
ensembl2hgnc.matrix <- function(df) {
  return(ensembl2hgnc.data.frame(df))
}

#' @rdname convert
#' @export
ensembl2hgnc.list <- function(df) {
  exp.list <- list()
  for (i in 1:length(df)) {
    exp.list[[i]] <- ensembl2hgnc(df[[i]])
  }
  return(exp.list)
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

#' Rename rownames of lists
#'
#' @description Rename colnames or rownames of list.
#' @export
#' @param x input list to rename element for.
#' @param ref.list List to take names from.
#' @param type Character specifying which names of x to change.
#' @param ref.list.type Character specifying which names of ref.list to change.
rename_list <- function(x, ref.list, type = "rownames", ref.list.type = "colnames") {
  renamed.list <- list()
  for (i in 1:length(x)) {
    df <- x[[i]]
    if (type == "rownames") {
      if (ref.list.type == "colnames") {
        rownames(df) <- colnames(ref.list[[i]])
      } else if (ref.list.type == "rownames") {
        rownames(df) <- rownames(ref.list[[i]])
      }
    } else if (type == "colnames") {
      if (ref.list.type == "colnames") {
        colnames(df) <- colnames(ref.list[[i]])
      } else if (ref.list.type == "rownames") {
        colnames(df) <- rownames(ref.list[[i]])
      }
    }
    renamed.lsit[[i]] <- df
  }
  return(renamed.list)
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
