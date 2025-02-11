// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// compute_forward_prob_cpp
NumericMatrix compute_forward_prob_cpp(const IntegerVector& x, const NumericVector& pInit, const NumericVector& pEmit, const NumericVector& Q, int p, int M, int K);
RcppExport SEXP _spacrt_compute_forward_prob_cpp(SEXP xSEXP, SEXP pInitSEXP, SEXP pEmitSEXP, SEXP QSEXP, SEXP pSEXP, SEXP MSEXP, SEXP KSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerVector& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type pInit(pInitSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type pEmit(pEmitSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type Q(QSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type M(MSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_forward_prob_cpp(x, pInit, pEmit, Q, p, M, K));
    return rcpp_result_gen;
END_RCPP
}
// compute_backward_prob_cpp
NumericMatrix compute_backward_prob_cpp(const IntegerVector& x, const NumericVector& pInit, const NumericVector& pEmit, const NumericVector& Q, int p, int M, int K);
RcppExport SEXP _spacrt_compute_backward_prob_cpp(SEXP xSEXP, SEXP pInitSEXP, SEXP pEmitSEXP, SEXP QSEXP, SEXP pSEXP, SEXP MSEXP, SEXP KSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerVector& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type pInit(pInitSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type pEmit(pEmitSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type Q(QSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type M(MSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_backward_prob_cpp(x, pInit, pEmit, Q, p, M, K));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_spacrt_compute_forward_prob_cpp", (DL_FUNC) &_spacrt_compute_forward_prob_cpp, 7},
    {"_spacrt_compute_backward_prob_cpp", (DL_FUNC) &_spacrt_compute_backward_prob_cpp, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_spacrt(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
