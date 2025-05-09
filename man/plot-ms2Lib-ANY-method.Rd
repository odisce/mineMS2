% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ms2Lib.R
\name{plot,ms2Lib,ANY-method}
\alias{plot,ms2Lib,ANY-method}
\title{Plot an element given an idx}
\usage{
\S4method{plot}{ms2Lib,ANY}(x, y, title = NULL, tkplot = FALSE, ...)
}
\arguments{
\item{x}{An ms2Lib oject.}

\item{y}{The index, a string starting by S if it a spectrum, P if it's a pattern or D if it's a dag.}

\item{title}{The title of the plot. If null, a title is automatically furnished.}

\item{tkplot}{if TRUE, for fragmentation graphs ("D"), the plot is dynamic and displayed using tkplot, otherwise it is static.}

\item{...}{supplementary arguments to be passed by the method.}
}
\value{
a fragPattern object or nothing
}
\description{
The method depends of the type of the index furnished. The following prefixes are supported:
\code{'P'} for a pattern, the plot method of fragPattern object is called (using the package igraph).
\code{'D'} for a DAG, the plot_dag method is called using the package igraph.
\code{'S'} for a spectrum, a spectrum is plotted calling the Plot method of Spectrum2 object.
Any other value will be removed.
}
\examples{
#' #Loading the data
data(m2l)

#Plotting a pattern
plot(m2l,"S10")

#The associated DAG
plot(m2l,"D10")

#Plotting a pattern
plot(m2l,"P53")
}
