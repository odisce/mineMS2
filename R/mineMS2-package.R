#' Mining of MS-MS spectra by FSM
#'
#' Mining of MS-MS spectra by frequent subgraph mining method.
#'
#' Mining of MS-MS spectra by frequent subgraph mining method followed. T
#'
#' @docType package
#' @importFrom MSnbase readMgfData precursorMz mz intensity
#' @importFrom igraph graph.empty set_vertex_attr add_edges adjacent_vertices V E get.edge.ids edge_attr
#' @import Rcpp
#' @name mineMS2-package
#' @useDynLib mineMS2
NULL