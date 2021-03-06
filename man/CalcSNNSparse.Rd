% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SNN.R
\name{CalcSNNSparse}
\alias{CalcSNNSparse}
\title{Function to convert the knn graph into the snn graph. Stored in a sparse
representation.}
\usage{
CalcSNNSparse(cell.names, k.param, nn.large, nn.ranked, prune.SNN,
  print.output)
}
\arguments{
\item{cell.names}{A vector of cell names which will correspond to the row/
column names of the SNN}

\item{k.param}{Defines nearest neighborhood when computing NN graph}

\item{nn.large}{Full KNN graph (computed with get.knn with k set to
k.param * k.scale)}

\item{nn.ranked}{Subset of full KNN graph that only contains the first
k.param nearest neighbors. Used to define Jaccard
distances between any two cells}

\item{prune.SNN}{Sets the cutoff for acceptable Jaccard distances when
computing the neighborhood overlap for the SNN
construction. Any edges with values less than or equal to
this will be set to 0 and removed from the SNN graph.
Essentially sets the strigency of pruning (0 --- no
pruning, 1 --- prune everything).}

\item{print.output}{Whether or not to print output to the console}
}
\value{
Returns an adjacency matrix representation of the SNN
                     graph
}
\description{
Function to convert the knn graph into the snn graph. Stored in a sparse
representation.
}
