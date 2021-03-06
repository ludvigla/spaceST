% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utility_functions.R
\name{convert}
\alias{convert}
\alias{ensembl2hgnc}
\alias{ensembl2hgnc.default}
\alias{ensembl2hgnc.list}
\alias{ensembl2hgnc.spaceST}
\title{Convert from ENSEMBL ID to HGNC symbol or MGI symbol}
\usage{
ensembl2hgnc(object, organism = "human")

\method{ensembl2hgnc}{default}(object, organism = "human")

\method{ensembl2hgnc}{list}(object, organism = "human")

\method{ensembl2hgnc}{spaceST}(object, organism = "human")
}
\arguments{
\item{object}{Data.frame, matrix or spaceST object with ENSEMBL ids as rownames.}

\item{organism}{Select organism database (human, mouse)}
}
\value{
Data.frame or matrix with converted names.
}
\description{
Function used to convert ENSEMBL gene ids of an expression matrix into HGNC/MGI symbols.
Duplicated gene names will be aggregated.
}
