% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/allClasses.R
\docType{class}
\name{ms2Lib-class}
\alias{ms2Lib-class}
\title{A class storing a set of MS/MS spectra, their associated graphs,
and some supplementary information.}
\description{
A class storing a set of MS/MS spectra, their associated graphs,
and some supplementary information.
}
\section{Slots}{

\describe{
\item{\code{spectraInfo}}{A data.frame including the spectra information at minima the "mz.precursor" field
and informations added by the user eventually.}

\item{\code{spectra}}{A list of spectra stored under the form of Spectrum2 object.}

\item{\code{ids}}{A chracter vector containing the ids of spectra.}

\item{\code{dags}}{A list storing the set of graphs object corresponding to the MS-MS spectra.}

\item{\code{losses}}{The discretized edge labels.}

\item{\code{fragments}}{The discretized nodes labels.}

\item{\code{patterns}}{A list storing the fragPattern objects.}

\item{\code{atoms}}{A list of the atoms used to build the formula as well as their maximum number.}

\item{\code{loss}}{A boolean indicating if the object is built with losses or fragments.}

\item{\code{reducedPatterns}}{The ids of the reduced patterns.}

\item{\code{k}}{The maximum depth of the constructed k path tree.}
}}

