#' Run DE analysis using the edgeR package to detect differentially expressed genes in each cluster
#'
#' @description Run DE analysis for each cluster of spots vs remaining spots
#' @export
#' @param spST spaceST object
#' @param clusters Integer vector of cluster IDs or character specifying existing cluster vector ("topic", "SNN")
#' @param verbose Print messages
#' @param pvalue pvalue cutoff
#' @importFrom scater newSCESet
#' @importFrom scran convertTo
#' @importFrom edgeR glmLRT glmFit topTags estimateDisp glmFit
#' @return Topic matrix
RunDEspaceST <- function(spST, clusters = "topic", verbose = TRUE, pvalue = 0.05) {
  stopifnot(length(spST@normalized) > 0)
  if (clusters == "topic") {
    stopifnot(!is.null(spST@meta.data$clusters))
    clusters <- spST@meta.data$clusters
    if (verbose) {
      message("Using clusters from hierarchical clustering of topic matrix.")
    }
  } else if (clusters == "SNN") {
    stopifnot(!is.null(spST@meta.data$clusters.snn))
    clusters <- spST@meta.data$clusters.snn
    if (verbose) {
      message("Using clusters from network based clustering of topic matrix.")
    }
  }
  cluster = factor(clusters)
  design.list <- list()
  tmp.cluster <- rep("NA", length(cluster))
  tmp.cluster.list <- list()
  for (i in 1:length(unique(cluster))){
    c <- sort(unique(cluster))[i]
    tmp.cluster[which(cluster %in% c)] <- 1
    tmp.cluster[which(cluster != c)] <- 2
    tmp.cluster.list[[i]] <- tmp.cluster
    design.list[[i]] <- model.matrix(~0 + tmp.cluster)
  }

  # This step will take some time to compute. Grab a coffee!
  fit.list <- parallel::mclapply(1:length(design.list), function(x) {
    design <- design.list[[x]]
    sce = newSCESet(countData = spST@normalized)
    y = convertTo(sce, type = "edgeR")
    y = estimateDisp(y, design)
    return(y)
  })

  # Intitiate empty lists to save reuslts in
  res.total.list <- parallel::mclapply(1:length(design.list), function(x) {
    # Obtain design matrix
    design <- design.list[[x]]
    # Obtain fitted model
    fit <- fit.list[[x]]
    # Name clusters
    clust.A = 1
    clust.B = 2
    contrast = numeric(ncol(design))
    contrast[clust.A] = 1
    contrast[clust.B] = -1
    res = glmLRT(fit, contrast = contrast)
    res.deg.final <- topTags(res, n = nrow(res$table), adjust.method = "BH", sort.by = "PValue", p.value = pvalue)$table
    #rownames(res.deg) <- rownames(df)
    if (is.null(res.deg.final)) {
      next
    }
    res.deg.final <- cbind(gene = rownames(res.deg.final), res.deg.final)
    colnames(res.deg.final) <- c("gene", colnames(res.deg.final))
    return(res.deg.final)
  })

  spST@DE.list <- res.total.list
  return(spST)
}


#' Volcano plot
#'
#' @description create a volcano plot of differentiually expressed genes
#' @param object spaceST object with differentially expressed genes
#' @param FDR.cutoff set cutoff threshold for adjusted p-values
#' @param logFC.cutoff set log2-foldchange cutoff
#' @param annotate.top.genes specify whether the top genes should be annotated in the plot
#' @param widget save html widget to file path
#' @importFrom plotly plot_ly layout as.widget
#' @importFrom htmlwidgets saveWidget
#' @export
volcano_plot_spaceST <- function(object, cluster, FDR.cutoff = 0.05, logFC.cutoff = 1, annotate.top.genes = NULL, widget.out = NULL) {
  stopifnot(class(object) == "spaceST",
            length(object@DE.list) > 0)
  DE.matrix <- object@DE.list[[cluster]]
  DE.matrix["group"] <- "NotSignificant"

  # change the grouping for the entries with significance but not a large enough Fold change
  DE.matrix[which(DE.matrix['FDR'] < FDR.cutoff & abs(DE.matrix['logFC']) > logFC.cutoff ),"group"] <- "Significant"

  if (!is.null(annotate.top.genes)) {
    top_peaks <- DE.matrix[with(DE.matrix, order(logFC, FDR)),][1:annotate.top.genes,]
    top_peaks <- rbind(top_peaks, DE.matrix[with(DE.matrix, order(-logFC, FDR)),][1:annotate.top.genes,])
  }

  a <- list()
  for (i in seq_len(nrow(top_peaks))) {
    m <- top_peaks[i, ]
    a[[i]] <- list(
      x = m[["logFC"]],
      y = -log10(m[["FDR"]]) + 0.5,
      text = m[["gene"]],
      xref = "x",
      yref = "y",
      showarrow = FALSE,
      arrowhead = 0,
      ax = 20,
      ay = -40
    )
  }

  # make the Plot.ly plot
  p <- plot_ly(data = DE.matrix,
               x = ~logFC,
               y = ~-log10(FDR),
               text = ~gene,
               mode = "markers",
               color = ~group, colors = c("grey", "red"),
               type = 'scatter') %>%
    layout(title ="Volcano Plot") %>%
    layout(annotations = a, xaxis = list(
      title = "log2(foldchange)"),
      yaxis = list(title = "-log10(adj. P-value)"))
  if (!is.null(widget.out)) {
    saveWidget(as.widget(p), widget.out)
  } else {
    p
  }
}


#' Plot DE heatmap
#'
#' @description Function used to plot heatmaps of DE genes expression
#' @export
#' @param object spaceST object with DE results
#' @param FDR.cutoff cutoff for adjusted p-values
#' @param logFC.cutoff cutof for logFC
#' @param log.2 log2 transform data
#' @param type select data type from spaceST object
#' @param scale.data logical speficying whether or not values should be scaled
#' @importFrom gplots heatmap.2
#' @importFrom dplyr group_by filter
#' @importFrom magrittr %>%
#' @return Heatmap of DE genes.
DE_heatmap <- function(object,
                       FDR.cutoff = 0.05,
                       logFC.cutoff = 2,
                       log.2 = T,
                       type = "normalized",
                       include.downregulated = F,
                       scale.data = T){
  stopifnot(class(object) == "spaceST",
            length(object@DE.list) > 0)
  if (type == "raw") {
    df <- object@expr
  } else if (type == "corrected") {
    df <- object@corrected
  } else if (type == "normalized") {
    df <- object@normalized
  }

  DE.genes <- do.call(rbind, object@DE.list)
  DE.genes$cluster <- rep(1:length(unique(object@meta.data$clusters)), each = dim(df)[1])
  DE.genes <- DE.genes %>% group_by(cluster) %>% dplyr::filter(FDR < FDR.cutoff)

  if (include.downregulated) {
    DE.genes <- DE.genes %>% dplyr::filter(abs(logFC) > logFC.cutoff)
  } else {
    DE.genes <- DE.genes %>% dplyr::filter(logFC > logFC.cutoff)
  }
  df <- df[as.character(DE.genes$gene), ]
  if (log.2) {
    df <- log2(df + 1)
    if (scale.data) {
      df <- scale(df, center = T, scale = T)
    }
  }

  # Order df by cluster
  clusters = object@meta.data$clusters
  ind <- order(sprintf("%03d", clusters))
  df <- df[, ind]
  clusters <- clusters[ind]

  mycol = colorRampPalette(c("dark blue", "cyan", "yellow", "red"))(256)
  #colors = c("#114477", "#4477AA", "#77AADD", "#117755", "#44AA88", "#99CCBB", "#777711", "#AAAA44", "#DDDD77", "#771111", "#AA4444", "#DD7777", "#771144", "#AA4477", "#DD77AA")
  colors = c("royalblue3",
             "mediumpurple4",
             "navajowhite2",
             "chocolate",
             "firebrick",
             "yellow2",
             "aquamarine",
             "orange1",
             "olivedrab2",
             "darkgreen",
             "pink",
             "black",
             "navy",
             "khaki3",
             "lightsteelblue1")
  clusters.col = colors[clusters]
  heatmap.2(x = as.matrix(df),
    key = TRUE,
    key.xlab = 'Identity',
    key.ylab = '',
    xlab = "feature",
    ylab = "gene expression",
    density.info = 'none',
    scale = "none",
    trace = 'none',
    symbreaks = F,
    revC = F,
    cexRow = 0.3,
    cexCol = 0.35,
    symkey = 0,
    dendrogram = 'row',
    hclustfun = function(m)hclust(m, method = 'ward.D2'),
    distfun = function(m)dist(m, method = 'euclidean'),
    col = mycol,
    Rowv = T,
    Colv = F,
    ColSideColors = clusters.col)
}


#' Detect variable genes using trendsceek package
#'
#' @description Find most variable genes using calc_varstats function from trendsceek.
#' @param x spaceST object
#' @param quantile.cutoff set cutoff level
#' @param method method used for linear regression
#' @importFrom trendsceek calc_varstats genefilter_exprmat
#' @export
var.genes <- function(
  x,
  quantile.cutoff = 0.9,
  method = 'glm',
  min.ncells.expr = 3,
  min.expr = 5
) {
  stopifnot(class(x) == "spaceST" & length(x@normalized) > 0)
  x <- x@normalized
  counts_filt = genefilter_exprmat(as.matrix(x), min.expr, min.ncells.expr)
  vargenes_stats = calc_varstats(counts_filt, counts_filt, quant.cutoff = quantile.cutoff, method = method)
  return(vargenes_stats)
}
