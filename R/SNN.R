#' Build SNN graph
#'
#' This function is a modified version of the BuildSNN function in the single cell genomics
#' tool \href{https://github.com/satijalab/seurat}{Seurat} (Version: 2.1.0, date: 2017-10-11).
#'
#' @param object Topic matrix obtained with compute.lda() from cellTree package
#' @param k.param Defines k for the k-nearest neighbor algorithm
#' @param k.scale Granularity option for k.param
#' @param plot.SNN Plot the SNN graph
#' @param prune.SNN Sets the cutoff for acceptable Jaccard distances when
#' computing the neighborhood overlap for the SNN construction. Any edges with
#' values less than or equal to this will be set to 0 and removed from the SNN
#' graph. Essentially sets the strigency of pruning (0 --- no pruning, 1 ---
#' prune everything).
#' @param print.output Whether or not to print output to the console
#' @param cluster.method Chose method to use for clustering. "topic" will run
#' topic_clusters() and cluster features based on topic porportion using hierarchical
#' clustering with euclidean distances and ward.D2. Clusters are chosen using the
#' unsupervised cutreeDynamic function from dynamicTreeCut package. "SNN" will run
#' a network graph based algorithm using methods from the Seurat package.
#' @param algorithm Algorithm for modularity optimization (1 = original Louvain
#' algorithm; 2 = Louvain algorithm with multilevel refinement; 3 = SLM algorithm)
#' @param perplexity Set perplexity for tsne if plot.SNN is TRUE
#' @param clusters Clusters used to color vertices in SNN ntwork graph
#' @param seed Set random seed for t-SNE algorithm
#' @param color.by.sample Color vertices by sample in network graph
#' @param tsne Two column matrix with tsne results, alternatively run tsne using preset parameters
#' @param layout Set custom layout of vertices in network graph. Will override t-SNE results
#' @param ... Parameters passed to RunModularityClusteringspaceST() function used for clustering
#' @importFrom FNN get.knn
#' @importFrom igraph plot.igraph categorical_pal
#' @importFrom graphics legend
#' @importFrom Matrix sparseMatrix
#' @return SNN matrix
#' @rdname BuildTopicSNN
#' @export
BuildTopicSNN <- function(
  object,
  k.param = 30,
  k.scale = 25,
  plot.SNN = FALSE,
  prune.SNN = 1/15,
  print.output = TRUE,
  clusters = NULL,
  seed = 0,
  cluster.method = "SNN",
  perplexity = 15,
  color.by.sample,
  tsne = NULL,
  layout = NULL,
  algorithm = 1,
  ...
) {
  if (!length(object@topics) > 0) {
    stop("Error: missing topics matrix. Run topic_compute to obtain topic matrix.")
  }
  data.use <- object@topics
  n.cells <- nrow(x = data.use)
  if (n.cells < k.param) {
    warning("k.param set larger than number of cells. Setting k.param to number of cells - 1.")
    k.param <- n.cells - 1
  }
  # find the k-nearest neighbors for each feature
  my.knn <- get.knn(
    data <- as.matrix(x = data.use),
    k = min(k.scale * k.param, n.cells - 1)
  )
  nn.ranked <- cbind(1:n.cells, my.knn$nn.index[, 1:(k.param-1)])
  nn.large <- my.knn$nn.index

  w <- CalcSNNSparse(
    cell.names = rownames(data.use),
    k.param = k.param,
    nn.large = nn.large,
    nn.ranked = nn.ranked,
    prune.SNN = prune.SNN,
    print.output = print.output
  )

  if (plot.SNN) {
      net <- graph.adjacency(
        adjmatrix = as.matrix(w),
        mode = "undirected",
        weighted = TRUE,
        diag = FALSE
      )
      if (is.null(tsne)) {
        set.seed(seed)
        tsne <- Rtsne(as.matrix(data.use),
                         dims = 2,
                         initial_dims = 50,
                         theta = 0.0,
                         check_duplicates = FALSE,
                         pca = TRUE,
                         perplexity = perplexity,
                         max_iter = 1000,
                         verbose = F)$Y
      }
      if (is.null(clusters)) {
        if (cluster.method == "topic") {
          clusters <- topic_clusters(
            omega = data.use
          )
        } else if (cluster.method == "SNN") {
          clusters <- RunModularityClusteringspaceST(
            SNN = w,
            print.output = print.output,
            set.algorithm = algorithm,
            ...
          )
          object@meta.data$clusters.snn <- clusters
        }
      } else if (is.null(clusters) & color.by.sample) {
        clusters <- object@coordinates$replicate
      }
      if (is.null(layout)) {
        layout = tsne
      }
      plot.igraph(
        x = net,
        layout = as.matrix(x = layout),
        edge.width = E(graph = net)$weight*5,
        vertex.label = NA,
        vertex.size = 3,
        vertex.color = clusters
      )
      legend(
        "topleft",
        legend = levels(factor(clusters)),
        pch = 16,
        col = categorical_pal(n = length(unique(clusters))),
        title = "cluster", cex = 0.75
      )
  }
  object@snn <- w
  return(object)
}

#' Runs the modularity optimizer java program (ModularityOptimizer.jar)
#'
#' @param SNN SNN              matrix to use as input for the clustering
#'                             algorithms
#' @param modularity           Modularity function to use in clustering (1 =
#'                             standard; 2 = alternative).
#' @param resolution           Value of the resolution parameter, use a value
#'                             above (below) 1.0 if you want to obtain a larger
#'                             (smaller) number of communities.
#' @param algorithm            Algorithm for modularity optimization (1 =
#'                             original Louvain algorithm; 2 = Louvain algorithm
#'                             with multilevel refinement; 3 = SLM algorithm)
#' @param n.start              Number of random starts.
#' @param n.iter               Maximal number of iterations per random start.
#' @param random.seed          Seed of the random number generator
#' @param print.output         Whether or not to print output to the console
#' @param temp.file.location   Directory where intermediate files will be written.
#' @return                     Seurat object with identities set to the results
#'                             of the clustering procedure.
#'
#' @importFrom utils read.table write.table
RunModularityClusteringspaceST <- function(
  SNN = NULL,
  set.modularity = 1,
  set.resolution = 0.8,
  set.algorithm = 1,
  n.start = 100,
  n.iter = 10,
  random.seed = 0,
  print.output = FALSE,
  temp.file.location = NULL
) {
  stopifnot(!is.null(SNN))
  seurat.dir <- system.file(package = "Seurat")
  ModularityJarFile <- paste0(seurat.dir, "/java/ModularityOptimizer.jar")
  seurat.dir.base <- strsplit(x = seurat.dir, split = "/")[[1]]
  seurat.dir <- paste0(
    seurat.dir.base[0:(length(x = seurat.dir.base) - 1)],
    collapse = "/"
  )
  seurat.dir <- paste0(seurat.dir, "/")
  diag(x = SNN) <- 0
  if (is.object(x = SNN)) {
    SNN <- as(object = SNN, Class = "dgTMatrix")
    edge <- cbind(i = SNN@j, j = SNN@i, x = SNN@x)
  }
  rownames(x = edge) <- NULL
  colnames(x = edge) <- NULL
  edge <- edge[! duplicated(x = edge[, 1:2]), ]
  temp.file.location <- SetIfNull(x = temp.file.location, default = seurat.dir)
  unique_ID <- sample(x = 10000:99999, size = 1)
  edge_file <- paste0(temp.file.location, "edge_", unique_ID, ".txt")
  output_file <- paste0(temp.file.location, "output_", unique_ID, ".txt")
  while (file.exists(edge_file)) {
    unique_ID <- sample(x = 10000:99999, size = 1)
    edge_file <- paste0(temp.file.location, "edge_", unique_ID, ".txt")
    output_file <- paste0(temp.file.location, "output", unique_ID, ".txt")
  }
  if (print.output) {
    print.output <- 1
  } else {
    print.output <- 0
  }
  write.table(
    x = edge,
    file = edge_file,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE
  )
  if (set.modularity == 2 && set.resolution > 1) {
    stop("error: resolution<1 for alternative modularity")
  }
  print(class(set.modularity))
  print(class(set.resolution))
  print(class(set.algorithm))
  print(class(n.start))
  print(class(n.iter))
  print(class(random.seed))
  print(class(print.output))
  command <- paste(
    "java -jar",
    shQuote(string = ModularityJarFile),
    shQuote(string = edge_file),
    shQuote(string = output_file),
    set.modularity,
    set.resolution,
    set.algorithm,
    n.start,
    n.iter,
    random.seed,
    print.output
  )
  system(command, wait = TRUE)
  ident.use <- read.table(file = output_file, header = FALSE, sep = "\t")[, 1]

  file.remove(edge_file)
  file.remove(output_file)
  return (ident.use + 1)
}

#' Helper function
#'
#' @param x R object
#' @param default R object to replace x with if x is NULL
#' @return default if x is NULL
#' @export
SetIfNull <- function (x, default)
{
  if (is.null(x = x)) {
    return(default)
  }
  else {
    return(x)
  }
}

#' Function to convert the knn graph into the snn graph. Stored in a sparse
#' representation.
#' @param cell.names    A vector of cell names which will correspond to the row/
#'                      column names of the SNN
#' @param k.param       Defines nearest neighborhood when computing NN graph
#' @param nn.large      Full KNN graph (computed with get.knn with k set to
#'                      k.param * k.scale)
#' @param nn.ranked     Subset of full KNN graph that only contains the first
#'                      k.param nearest neighbors. Used to define Jaccard
#'                      distances between any two cells
#' @param prune.snn     Sets the cutoff for acceptable Jaccard distances when
#'                      computing the neighborhood overlap for the SNN
#'                      construction. Any edges with values less than or equal to
#'                      this will be set to 0 and removed from the SNN graph.
#'                      Essentially sets the strigency of pruning (0 --- no
#'                      pruning, 1 --- prune everything).
#' @param print.output  Whether or not to print output to the console
#' @return              Returns an adjacency matrix representation of the SNN
#'                      graph
#'
#' @importFrom utils txtProgressBar setTxtProgressBar
CalcSNNSparse <- function (cell.names, k.param, nn.large, nn.ranked, prune.SNN,
          print.output)
{
  n.cells <- length(cell.names)
  counter <- 1
  idx1 <- vector(mode = "integer", length = n.cells^2/k.param)
  idx2 <- vector(mode = "integer", length = n.cells^2/k.param)
  edge.weight <- vector(mode = "double", length = n.cells^2/k.param)
  id <- 1
  if (print.output) {
    print("Constructing SNN")
    pb <- txtProgressBar(min = 0, max = n.cells, style = 3)
  }
  for (i in 1:n.cells) {
    for (j in 1:ncol(x = nn.large)) {
      s <- intersect(x = nn.ranked[i, ], y = nn.ranked[nn.large[i,
                                                                j], ])
      u <- union(nn.ranked[i, ], nn.ranked[nn.large[i,
                                                    j], ])
      e <- length(x = s)/length(x = u)
      if (e > prune.SNN) {
        idx1[id] <- i
        idx2[id] <- nn.large[i, j]
        edge.weight[id] <- e
        id <- id + 1
      }
    }
    if (print.output) {
      setTxtProgressBar(pb = pb, value = i)
    }
  }
  if (print.output) {
    close(con = pb)
  }
  idx1 <- idx1[!is.na(x = idx1) & idx1 != 0]
  idx2 <- idx2[!is.na(x = idx2) & idx2 != 0]
  edge.weight <- edge.weight[!is.na(x = edge.weight) & edge.weight !=
                               0]
  w <- sparseMatrix(i = idx1, j = idx2, x = edge.weight, dims = c(n.cells,
                                                                  n.cells))
  diag(x = w) <- 1
  rownames(x = w) <- cell.names
  colnames(x = w) <- cell.names
  return(w)
}
