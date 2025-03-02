% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lbvem.R
\name{lbvem}
\alias{lbvem}
\title{Co-clustering main function using gaussian LBVEM}
\usage{
lbvem(data, nbcoclust, init = 10)
}
\arguments{
\item{data}{continuous data (matrix or dataframe)}

\item{nbcoclust}{2d vector specifying the number of clusters (row and column respectively)}

\item{init}{number of iterations for the initialization}
}
\value{
a coclust result object
}
\description{
Co-clustering main function using gaussian LBVEM
}
\examples{
library(coclust)
library(blockcluster)
data("gaussiandata")
lbvem(gaussiandata, c(3, 3), init=50)
}
