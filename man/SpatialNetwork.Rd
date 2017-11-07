% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SpatialNetwork.R
\name{SpatialNetwork}
\alias{SpatialNetwork}
\title{spaceST network graphs}
\usage{
SpatialNetwork(object, subset_by = NULL, color_by = "clusters",
  select.clust.group = 1, select.rep.group = NULL,
  mode = "fruchtermanreingold", vertex.size = 1, edge.color = "black",
  edge.size = 0.1, clusters.snn = F, ...)
}
\arguments{
\item{subset_by}{Select attribute to subset igraph by, default NULL (options: "clusters", "replicate")}

\item{color_by}{Select attribute to color vertices by (options: "clusters", "replicate", "clust_and_rep)}

\item{select.clust.group}{Select cluster by cluster number}

\item{select.rep.group}{Select replicate by replicate ID}

\item{mode}{Select mode for network object, default NULL (options: "spatial", "circle", "kamadakawai", ...
For more options, see available placement algorithms in the sna package)}

\item{...}{Parameters passed to ggnet2}
}
\description{
Visualize SNN networks for spatial data
}