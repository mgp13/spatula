// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// triplets_to_pairs
IntegerMatrix triplets_to_pairs(const IntegerMatrix& triplets);
RcppExport SEXP _gagarin_triplets_to_pairs(SEXP tripletsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerMatrix& >::type triplets(tripletsSEXP);
    rcpp_result_gen = Rcpp::wrap(triplets_to_pairs(triplets));
    return rcpp_result_gen;
END_RCPP
}
// dgCMat_to_pairs_cpp
IntegerMatrix dgCMat_to_pairs_cpp(const IntegerVector& ivec, const IntegerVector& pvec, const IntegerVector& xvec);
RcppExport SEXP _gagarin_dgCMat_to_pairs_cpp(SEXP ivecSEXP, SEXP pvecSEXP, SEXP xvecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerVector& >::type ivec(ivecSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type pvec(pvecSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type xvec(xvecSEXP);
    rcpp_result_gen = Rcpp::wrap(dgCMat_to_pairs_cpp(ivec, pvec, xvec));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_gagarin_triplets_to_pairs", (DL_FUNC) &_gagarin_triplets_to_pairs, 1},
    {"_gagarin_dgCMat_to_pairs_cpp", (DL_FUNC) &_gagarin_dgCMat_to_pairs_cpp, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_gagarin(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
