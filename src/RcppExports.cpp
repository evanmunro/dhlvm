// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// posteriorLikelihood
double posteriorLikelihood(NumericMatrix data, NumericVector groups, NumericMatrix pi, List beta);
RcppExport SEXP _dhlvm_posteriorLikelihood(SEXP dataSEXP, SEXP groupsSEXP, SEXP piSEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type data(dataSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type groups(groupsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type pi(piSEXP);
    Rcpp::traits::input_parameter< List >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(posteriorLikelihood(data, groups, pi, beta));
    return rcpp_result_gen;
END_RCPP
}
// hlc_cpp_noZ
List hlc_cpp_noZ(NumericMatrix data, NumericVector groups, List eta, NumericMatrix alpha, int steps);
RcppExport SEXP _dhlvm_hlc_cpp_noZ(SEXP dataSEXP, SEXP groupsSEXP, SEXP etaSEXP, SEXP alphaSEXP, SEXP stepsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type data(dataSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type groups(groupsSEXP);
    Rcpp::traits::input_parameter< List >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< int >::type steps(stepsSEXP);
    rcpp_result_gen = Rcpp::wrap(hlc_cpp_noZ(data, groups, eta, alpha, steps));
    return rcpp_result_gen;
END_RCPP
}
// hlc_cpp
List hlc_cpp(NumericMatrix data, NumericVector groups, List eta, NumericMatrix alpha, int steps);
RcppExport SEXP _dhlvm_hlc_cpp(SEXP dataSEXP, SEXP groupsSEXP, SEXP etaSEXP, SEXP alphaSEXP, SEXP stepsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type data(dataSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type groups(groupsSEXP);
    Rcpp::traits::input_parameter< List >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< int >::type steps(stepsSEXP);
    rcpp_result_gen = Rcpp::wrap(hlc_cpp(data, groups, eta, alpha, steps));
    return rcpp_result_gen;
END_RCPP
}
// dhlc_cpp_noZ
List dhlc_cpp_noZ(NumericMatrix data, NumericVector groups, List eta, int v0, int s0, double tune, int K, int T, int steps);
RcppExport SEXP _dhlvm_dhlc_cpp_noZ(SEXP dataSEXP, SEXP groupsSEXP, SEXP etaSEXP, SEXP v0SEXP, SEXP s0SEXP, SEXP tuneSEXP, SEXP KSEXP, SEXP TSEXP, SEXP stepsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type data(dataSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type groups(groupsSEXP);
    Rcpp::traits::input_parameter< List >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< int >::type v0(v0SEXP);
    Rcpp::traits::input_parameter< int >::type s0(s0SEXP);
    Rcpp::traits::input_parameter< double >::type tune(tuneSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type T(TSEXP);
    Rcpp::traits::input_parameter< int >::type steps(stepsSEXP);
    rcpp_result_gen = Rcpp::wrap(dhlc_cpp_noZ(data, groups, eta, v0, s0, tune, K, T, steps));
    return rcpp_result_gen;
END_RCPP
}
// dhlc_cpp
List dhlc_cpp(NumericMatrix data, NumericVector groups, List eta, int v0, int s0, double tune, int K, int T, int steps);
RcppExport SEXP _dhlvm_dhlc_cpp(SEXP dataSEXP, SEXP groupsSEXP, SEXP etaSEXP, SEXP v0SEXP, SEXP s0SEXP, SEXP tuneSEXP, SEXP KSEXP, SEXP TSEXP, SEXP stepsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type data(dataSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type groups(groupsSEXP);
    Rcpp::traits::input_parameter< List >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< int >::type v0(v0SEXP);
    Rcpp::traits::input_parameter< int >::type s0(s0SEXP);
    Rcpp::traits::input_parameter< double >::type tune(tuneSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type T(TSEXP);
    Rcpp::traits::input_parameter< int >::type steps(stepsSEXP);
    rcpp_result_gen = Rcpp::wrap(dhlc_cpp(data, groups, eta, v0, s0, tune, K, T, steps));
    return rcpp_result_gen;
END_RCPP
}
// sampleZ_old
List sampleZ_old(List data, NumericMatrix theta, NumericMatrix beta);
RcppExport SEXP _dhlvm_sampleZ_old(SEXP dataSEXP, SEXP thetaSEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type data(dataSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(sampleZ_old(data, theta, beta));
    return rcpp_result_gen;
END_RCPP
}
// sampleGamma_old
NumericMatrix sampleGamma_old(List z, NumericMatrix sigma, NumericMatrix gammaprev, double shrink);
RcppExport SEXP _dhlvm_sampleGamma_old(SEXP zSEXP, SEXP sigmaSEXP, SEXP gammaprevSEXP, SEXP shrinkSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type z(zSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type gammaprev(gammaprevSEXP);
    Rcpp::traits::input_parameter< double >::type shrink(shrinkSEXP);
    rcpp_result_gen = Rcpp::wrap(sampleGamma_old(z, sigma, gammaprev, shrink));
    return rcpp_result_gen;
END_RCPP
}
// sampleBetaLDA
NumericMatrix sampleBetaLDA(List data, List z, NumericMatrix eta);
RcppExport SEXP _dhlvm_sampleBetaLDA(SEXP dataSEXP, SEXP zSEXP, SEXP etaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type data(dataSEXP);
    Rcpp::traits::input_parameter< List >::type z(zSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type eta(etaSEXP);
    rcpp_result_gen = Rcpp::wrap(sampleBetaLDA(data, z, eta));
    return rcpp_result_gen;
END_RCPP
}
// discreteLDS_cpp
List discreteLDS_cpp(List data, NumericMatrix eta, double v0, double s0, double tune, int K, int V, int steps);
RcppExport SEXP _dhlvm_discreteLDS_cpp(SEXP dataSEXP, SEXP etaSEXP, SEXP v0SEXP, SEXP s0SEXP, SEXP tuneSEXP, SEXP KSEXP, SEXP VSEXP, SEXP stepsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type data(dataSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< double >::type v0(v0SEXP);
    Rcpp::traits::input_parameter< double >::type s0(s0SEXP);
    Rcpp::traits::input_parameter< double >::type tune(tuneSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type V(VSEXP);
    Rcpp::traits::input_parameter< int >::type steps(stepsSEXP);
    rcpp_result_gen = Rcpp::wrap(discreteLDS_cpp(data, eta, v0, s0, tune, K, V, steps));
    return rcpp_result_gen;
END_RCPP
}
// discreteMS_cpp
List discreteMS_cpp(NumericMatrix data, NumericMatrix eta, NumericMatrix alpha, int K, int steps);
RcppExport SEXP _dhlvm_discreteMS_cpp(SEXP dataSEXP, SEXP etaSEXP, SEXP alphaSEXP, SEXP KSEXP, SEXP stepsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type data(dataSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type steps(stepsSEXP);
    rcpp_result_gen = Rcpp::wrap(discreteMS_cpp(data, eta, alpha, K, steps));
    return rcpp_result_gen;
END_RCPP
}
// testRgamma
double testRgamma(double v0, double s0);
RcppExport SEXP _dhlvm_testRgamma(SEXP v0SEXP, SEXP s0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type v0(v0SEXP);
    Rcpp::traits::input_parameter< double >::type s0(s0SEXP);
    rcpp_result_gen = Rcpp::wrap(testRgamma(v0, s0));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_dhlvm_posteriorLikelihood", (DL_FUNC) &_dhlvm_posteriorLikelihood, 4},
    {"_dhlvm_hlc_cpp_noZ", (DL_FUNC) &_dhlvm_hlc_cpp_noZ, 5},
    {"_dhlvm_hlc_cpp", (DL_FUNC) &_dhlvm_hlc_cpp, 5},
    {"_dhlvm_dhlc_cpp_noZ", (DL_FUNC) &_dhlvm_dhlc_cpp_noZ, 9},
    {"_dhlvm_dhlc_cpp", (DL_FUNC) &_dhlvm_dhlc_cpp, 9},
    {"_dhlvm_sampleZ_old", (DL_FUNC) &_dhlvm_sampleZ_old, 3},
    {"_dhlvm_sampleGamma_old", (DL_FUNC) &_dhlvm_sampleGamma_old, 4},
    {"_dhlvm_sampleBetaLDA", (DL_FUNC) &_dhlvm_sampleBetaLDA, 3},
    {"_dhlvm_discreteLDS_cpp", (DL_FUNC) &_dhlvm_discreteLDS_cpp, 8},
    {"_dhlvm_discreteMS_cpp", (DL_FUNC) &_dhlvm_discreteMS_cpp, 5},
    {"_dhlvm_testRgamma", (DL_FUNC) &_dhlvm_testRgamma, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_dhlvm(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
