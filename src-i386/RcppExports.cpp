// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// closeMatch
IntegerVector closeMatch(NumericVector x, NumericVector y, IntegerVector xidx, IntegerVector yidx, double ppm, double dmz);
RcppExport SEXP _mineMS2_closeMatch(SEXP xSEXP, SEXP ySEXP, SEXP xidxSEXP, SEXP yidxSEXP, SEXP ppmSEXP, SEXP dmzSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type xidx(xidxSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type yidx(yidxSEXP);
    Rcpp::traits::input_parameter< double >::type ppm(ppmSEXP);
    Rcpp::traits::input_parameter< double >::type dmz(dmzSEXP);
    rcpp_result_gen = Rcpp::wrap(closeMatch(x, y, xidx, yidx, ppm, dmz));
    return rcpp_result_gen;
END_RCPP
}
// decomposeMass
List decomposeMass(double mass, double tolerance, int b, std::vector<double> components);
RcppExport SEXP _mineMS2_decomposeMass(SEXP massSEXP, SEXP toleranceSEXP, SEXP bSEXP, SEXP componentsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type mass(massSEXP);
    Rcpp::traits::input_parameter< double >::type tolerance(toleranceSEXP);
    Rcpp::traits::input_parameter< int >::type b(bSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type components(componentsSEXP);
    rcpp_result_gen = Rcpp::wrap(decomposeMass(mass, tolerance, b, components));
    return rcpp_result_gen;
END_RCPP
}
// formulaExtension
List formulaExtension(NumericVector masses, NumericVector mzlim, IntegerMatrix formula, LogicalVector hatoms, IntegerVector hhatom);
RcppExport SEXP _mineMS2_formulaExtension(SEXP massesSEXP, SEXP mzlimSEXP, SEXP formulaSEXP, SEXP hatomsSEXP, SEXP hhatomSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type masses(massesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mzlim(mzlimSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type formula(formulaSEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type hatoms(hatomsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type hhatom(hhatomSEXP);
    rcpp_result_gen = Rcpp::wrap(formulaExtension(masses, mzlim, formula, hatoms, hhatom));
    return rcpp_result_gen;
END_RCPP
}
// mineClosedDags
Rcpp::List mineClosedDags(List& vertices_list, List& edges_list, LogicalVector& processing, IntegerVector num, IntegerVector k, IntegerVector size_min, LogicalVector prec_only);
RcppExport SEXP _mineMS2_mineClosedDags(SEXP vertices_listSEXP, SEXP edges_listSEXP, SEXP processingSEXP, SEXP numSEXP, SEXP kSEXP, SEXP size_minSEXP, SEXP prec_onlySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List& >::type vertices_list(vertices_listSEXP);
    Rcpp::traits::input_parameter< List& >::type edges_list(edges_listSEXP);
    Rcpp::traits::input_parameter< LogicalVector& >::type processing(processingSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type num(numSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type k(kSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type size_min(size_minSEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type prec_only(prec_onlySEXP);
    rcpp_result_gen = Rcpp::wrap(mineClosedDags(vertices_list, edges_list, processing, num, k, size_min, prec_only));
    return rcpp_result_gen;
END_RCPP
}
// MST
DataFrame MST(DataFrame edges, NumericVector scores, int num_vertices);
RcppExport SEXP _mineMS2_MST(SEXP edgesSEXP, SEXP scoresSEXP, SEXP num_verticesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame >::type edges(edgesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type scores(scoresSEXP);
    Rcpp::traits::input_parameter< int >::type num_vertices(num_verticesSEXP);
    rcpp_result_gen = Rcpp::wrap(MST(edges, scores, num_vertices));
    return rcpp_result_gen;
END_RCPP
}
// scorePattern
double scorePattern(DataFrame edges, NumericVector scores, int num_vertices);
RcppExport SEXP _mineMS2_scorePattern(SEXP edgesSEXP, SEXP scoresSEXP, SEXP num_verticesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame >::type edges(edgesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type scores(scoresSEXP);
    Rcpp::traits::input_parameter< int >::type num_vertices(num_verticesSEXP);
    rcpp_result_gen = Rcpp::wrap(scorePattern(edges, scores, num_vertices));
    return rcpp_result_gen;
END_RCPP
}
// select_patterns_from_spectra
IntegerVector select_patterns_from_spectra(List pat_list, int sid);
RcppExport SEXP _mineMS2_select_patterns_from_spectra(SEXP pat_listSEXP, SEXP sidSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type pat_list(pat_listSEXP);
    Rcpp::traits::input_parameter< int >::type sid(sidSEXP);
    rcpp_result_gen = Rcpp::wrap(select_patterns_from_spectra(pat_list, sid));
    return rcpp_result_gen;
END_RCPP
}
// patterns_from_spectra
List patterns_from_spectra(List pat_list, int num_spectra);
RcppExport SEXP _mineMS2_patterns_from_spectra(SEXP pat_listSEXP, SEXP num_spectraSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type pat_list(pat_listSEXP);
    Rcpp::traits::input_parameter< int >::type num_spectra(num_spectraSEXP);
    rcpp_result_gen = Rcpp::wrap(patterns_from_spectra(pat_list, num_spectra));
    return rcpp_result_gen;
END_RCPP
}
// FindEqualGreaterM
IntegerVector FindEqualGreaterM(NumericVector inv, NumericVector values);
RcppExport SEXP _mineMS2_FindEqualGreaterM(SEXP invSEXP, SEXP valuesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type inv(invSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type values(valuesSEXP);
    rcpp_result_gen = Rcpp::wrap(FindEqualGreaterM(inv, values));
    return rcpp_result_gen;
END_RCPP
}
// findLimDensity
IntegerVector findLimDensity(NumericVector seq, int istart, int state);
RcppExport SEXP _mineMS2_findLimDensity(SEXP seqSEXP, SEXP istartSEXP, SEXP stateSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type seq(seqSEXP);
    Rcpp::traits::input_parameter< int >::type istart(istartSEXP);
    Rcpp::traits::input_parameter< int >::type state(stateSEXP);
    rcpp_result_gen = Rcpp::wrap(findLimDensity(seq, istart, state));
    return rcpp_result_gen;
END_RCPP
}
// formulaFromString
List formulaFromString(std::string formula, std::vector<std::string> names_atoms);
RcppExport SEXP _mineMS2_formulaFromString(SEXP formulaSEXP, SEXP names_atomsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type formula(formulaSEXP);
    Rcpp::traits::input_parameter< std::vector<std::string> >::type names_atoms(names_atomsSEXP);
    rcpp_result_gen = Rcpp::wrap(formulaFromString(formula, names_atoms));
    return rcpp_result_gen;
END_RCPP
}
// disjointBins
NumericVector disjointBins(NumericVector points, NumericVector lower_lim, NumericVector upper_lim, NumericVector mean_bin);
RcppExport SEXP _mineMS2_disjointBins(SEXP pointsSEXP, SEXP lower_limSEXP, SEXP upper_limSEXP, SEXP mean_binSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type points(pointsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lower_lim(lower_limSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type upper_lim(upper_limSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mean_bin(mean_binSEXP);
    rcpp_result_gen = Rcpp::wrap(disjointBins(points, lower_lim, upper_lim, mean_bin));
    return rcpp_result_gen;
END_RCPP
}
// checkInter
List checkInter(NumericVector a_min, NumericVector a_max, NumericVector b_min, NumericVector b_max);
RcppExport SEXP _mineMS2_checkInter(SEXP a_minSEXP, SEXP a_maxSEXP, SEXP b_minSEXP, SEXP b_maxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type a_min(a_minSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type a_max(a_maxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type b_min(b_minSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type b_max(b_maxSEXP);
    rcpp_result_gen = Rcpp::wrap(checkInter(a_min, a_max, b_min, b_max));
    return rcpp_result_gen;
END_RCPP
}
// find_combinations_ranges
List find_combinations_ranges(NumericVector bmin, NumericVector bmax, NumericVector cmzmax);
RcppExport SEXP _mineMS2_find_combinations_ranges(SEXP bminSEXP, SEXP bmaxSEXP, SEXP cmzmaxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type bmin(bminSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type bmax(bmaxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type cmzmax(cmzmaxSEXP);
    rcpp_result_gen = Rcpp::wrap(find_combinations_ranges(bmin, bmax, cmzmax));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_mineMS2_closeMatch", (DL_FUNC) &_mineMS2_closeMatch, 6},
    {"_mineMS2_decomposeMass", (DL_FUNC) &_mineMS2_decomposeMass, 4},
    {"_mineMS2_formulaExtension", (DL_FUNC) &_mineMS2_formulaExtension, 5},
    {"_mineMS2_mineClosedDags", (DL_FUNC) &_mineMS2_mineClosedDags, 7},
    {"_mineMS2_MST", (DL_FUNC) &_mineMS2_MST, 3},
    {"_mineMS2_scorePattern", (DL_FUNC) &_mineMS2_scorePattern, 3},
    {"_mineMS2_select_patterns_from_spectra", (DL_FUNC) &_mineMS2_select_patterns_from_spectra, 2},
    {"_mineMS2_patterns_from_spectra", (DL_FUNC) &_mineMS2_patterns_from_spectra, 2},
    {"_mineMS2_FindEqualGreaterM", (DL_FUNC) &_mineMS2_FindEqualGreaterM, 2},
    {"_mineMS2_findLimDensity", (DL_FUNC) &_mineMS2_findLimDensity, 3},
    {"_mineMS2_formulaFromString", (DL_FUNC) &_mineMS2_formulaFromString, 2},
    {"_mineMS2_disjointBins", (DL_FUNC) &_mineMS2_disjointBins, 4},
    {"_mineMS2_checkInter", (DL_FUNC) &_mineMS2_checkInter, 4},
    {"_mineMS2_find_combinations_ranges", (DL_FUNC) &_mineMS2_find_combinations_ranges, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_mineMS2(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}