# vi: fdm=marker ts=4 et cc=80 tw=80

# Constants {{{1
################################################################

TEST.DIR <- file.path(getwd(), '..')
RES.DIR  <- file.path(TEST.DIR, 'res')
OUT.DIR  <- file.path(TEST.DIR, 'output')
if ( ! dir.exists(OUT.DIR))
    dir.create(OUT.DIR)

# Test basic processing {{{1
################################################################################

test_basic_processing <- function() {
    
    path_mgf <- file.path(RES.DIR, 'basic.mgf')
    path_input_graph <- file.path(RES.DIR, 'basic.graphml')
    path_expected_output_graph <- file.path(RES.DIR, 'basic_annotated.graphml')
    path_output_graph <- file.path(OUT.DIR, 'test_basic_processing.graphml')

    m2l <- mineMS2::ms2Lib(path_mgf)
    testthat::expect_is(m2l, 'ms2Lib')
    infos <- mineMS2::getInfo(m2l, "S")
    ids <- paste(paste0("MZ", infos[,"mz.precursor"]), sep = "_")
    m2l <- mineMS2::setIds(m2l, ids)
    m2l <- mineMS2::discretizeMassLosses(m2l,
                                         dmz = 0.008,
                                         ppm = 8,
                                         heteroAtoms = FALSE,
                                         maxFrags = 15)
    m2l <- mineMS2::mineClosedSubgraphs(m2l,
                                        sizeMin = 1,
                                        count = 2)
    mineMS2::plotPatterns(m2l, full = TRUE)
    
    net_gnps <- igraph::read_graph(path_input_graph, "graphml")
    net_gnps <- igraph::simplify(net_gnps,
                                 remove.multiple = FALSE,
                                 edge.attr.comb = "ignore")
    testthat::expect_is(net_gnps, 'igraph')
    
    patterns <- mineMS2::findPatternsExplainingComponents(m2l,
                                                          components,
                                                          metric=c("recall",
                                                                   "precision",
                                                                   "size"))
    testthat::expect_is(patterns, 'list')
    
    annotated_net <- mineMS2::annotateNetwork(components,
                                              net_gnps,
                                              patterns)
    testthat::expect_is(annotated_net, 'igraph')

    igraph::write_graph(graph = annotated_net,
                        format = "graphml",
                        file = path_output_graph)

    # Compare outputs
    print(tools::md5sum(path_output_graph))
    print(tools::md5sum(path_expected_output_graph))
    testthat::expect_equal(tools::md5sum(path_output_graph),
                           tools::md5sum(path_expected_output_graph))
}

# Test MGF filename {{{1
################################################################################

test_basic_mgf_filename <- function() {
    
    path_mgf <- file.path(RES.DIR, 'basic.mgf')
    mgf_with_no_ext_mgf <- file.path(OUT.DIR, 'basic.dat')
    file.copy(path_mgf, mgf_with_no_ext_mgf)
    
    m2l <- mineMS2::ms2Lib(mgf_with_no_ext_mgf)
    testthat::expect_is(m2l, 'ms2Lib')
}

# Main {{{1
################################################################################

testthat::context('Main tests.')
testthat::test_that('A basic processing works fine.', test_basic_processing())
testthat::test_that('We can provide any filename for a single MGF file', test_basic_mgf_filename())
