% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spatial_tools.R
\name{spatial.heatmap}
\alias{spatial.heatmap}
\title{Visualize spatial heatmaps}
\usage{
spatial.heatmap(object, value, type = "expr", HE.list = NULL,
  palette = "spectral", invert.heatmap = FALSE, hide.legend = FALSE,
  bg.black = FALSE, xlim = c(1, 33), ylim = c(1, 35), ncols = NULL,
  arrange = TRUE, col.title = NULL, scale = "", size = 3,
  hide.dropouts = F, ...)
}
\arguments{
\item{object}{Object of class spaceST.}

\item{value}{Target vector to visualize. This value can be chosen from any slot of the spaceST object
containing data linked to array spots.}

\item{type}{Select dataset where the value can be found [options: "expr", "norm.data", "pca"]}

\item{HE.list}{List of paths to HE images in jpeg format that should be used as a background for the
spatial heatmap.}

\item{palette}{Color palette used for spatial heamtap [options: ""expr", "green.to.blue", "spectral", "offwhite.to.black", "
"viridis", "magma", "plasma", "cividis"]}

\item{invert.heatmap}{Invert color gradient.}

\item{hide.legend}{Exclude legend.}

\item{bg.black}{Set background color to black.}

\item{xlim, ylim}{Set limits of x/y axes. [default: xlim = c(1, 34), ylim = c(1, 36)]}

\item{ncols}{Number of columns in arranged plot table.}

\item{arrange}{Arrange plots.}

\item{col.title}{Give the color legend a title.}

\item{scale}{Set this parameter to "colwise" to scale each column separately.}

\item{size}{Set size of markers.}

\item{hide.dropouts}{Logical; Set to TRUE if you want to hide spots with 0 values.}

\item{...}{Parameters passed to geom_point.}
}
\description{
Generates a scatter plot of array spots colored by some value.
}
