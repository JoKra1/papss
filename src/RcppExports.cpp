// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// wrapAmSolve
Rcpp::List wrapAmSolve(const Eigen::Map<Eigen::MatrixXd> X, const Eigen::Map<Eigen::VectorXd> y, const Eigen::Map<Eigen::VectorXd>& initCf, const Rcpp::StringVector& constraints, const Rcpp::IntegerVector& lambdaTermFreq, const Rcpp::IntegerVector& lambdaTermDim, int startIndex, int maxIter, int maxIterOptim, double tol, bool shouldCollectProgress, double startLambda, bool shouldAccumulH);
RcppExport SEXP _papss_wrapAmSolve(SEXP XSEXP, SEXP ySEXP, SEXP initCfSEXP, SEXP constraintsSEXP, SEXP lambdaTermFreqSEXP, SEXP lambdaTermDimSEXP, SEXP startIndexSEXP, SEXP maxIterSEXP, SEXP maxIterOptimSEXP, SEXP tolSEXP, SEXP shouldCollectProgressSEXP, SEXP startLambdaSEXP, SEXP shouldAccumulHSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type X(XSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd> >::type y(ySEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd>& >::type initCf(initCfSEXP);
    Rcpp::traits::input_parameter< const Rcpp::StringVector& >::type constraints(constraintsSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type lambdaTermFreq(lambdaTermFreqSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type lambdaTermDim(lambdaTermDimSEXP);
    Rcpp::traits::input_parameter< int >::type startIndex(startIndexSEXP);
    Rcpp::traits::input_parameter< int >::type maxIter(maxIterSEXP);
    Rcpp::traits::input_parameter< int >::type maxIterOptim(maxIterOptimSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< bool >::type shouldCollectProgress(shouldCollectProgressSEXP);
    Rcpp::traits::input_parameter< double >::type startLambda(startLambdaSEXP);
    Rcpp::traits::input_parameter< bool >::type shouldAccumulH(shouldAccumulHSEXP);
    rcpp_result_gen = Rcpp::wrap(wrapAmSolve(X, y, initCf, constraints, lambdaTermFreq, lambdaTermDim, startIndex, maxIter, maxIterOptim, tol, shouldCollectProgress, startLambda, shouldAccumulH));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_papss_wrapAmSolve", (DL_FUNC) &_papss_wrapAmSolve, 13},
    {NULL, NULL, 0}
};

RcppExport void R_init_papss(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
