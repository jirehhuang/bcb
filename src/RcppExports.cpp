// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// obj_cpp
double obj_cpp(arma::vec par, arma::uvec pt_dim, arma::mat x_tp_, arma::vec b_t, arma::mat a_p);
RcppExport SEXP _bcb_obj_cpp(SEXP parSEXP, SEXP pt_dimSEXP, SEXP x_tp_SEXP, SEXP b_tSEXP, SEXP a_pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type par(parSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type pt_dim(pt_dimSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x_tp_(x_tp_SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type b_t(b_tSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type a_p(a_pSEXP);
    rcpp_result_gen = Rcpp::wrap(obj_cpp(par, pt_dim, x_tp_, b_t, a_p));
    return rcpp_result_gen;
END_RCPP
}
// con_cpp
Rcpp::List con_cpp(Rcpp::NumericVector par, arma::uvec pt_dim, arma::mat x_tp_, arma::vec b_t);
RcppExport SEXP _bcb_con_cpp(SEXP parSEXP, SEXP pt_dimSEXP, SEXP x_tp_SEXP, SEXP b_tSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type par(parSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type pt_dim(pt_dimSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x_tp_(x_tp_SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type b_t(b_tSEXP);
    rcpp_result_gen = Rcpp::wrap(con_cpp(par, pt_dim, x_tp_, b_t));
    return rcpp_result_gen;
END_RCPP
}
// beta_cpp
void beta_cpp(arma::mat& X, arma::vec& y, arma::vec& beta);
RcppExport SEXP _bcb_beta_cpp(SEXP XSEXP, SEXP ySEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type beta(betaSEXP);
    beta_cpp(X, y, beta);
    return R_NilValue;
END_RCPP
}
// lm_cpp
void lm_cpp(arma::mat& X, arma::vec& y, arma::vec& values);
RcppExport SEXP _bcb_lm_cpp(SEXP XSEXP, SEXP ySEXP, SEXP valuesSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type values(valuesSEXP);
    lm_cpp(X, y, values);
    return R_NilValue;
END_RCPP
}
// lm_nig
void lm_nig(arma::vec& Xty, arma::vec& m_0, arma::mat& Lambda_0, arma::mat& Lambda_n, arma::mat& V_n, double yty, double a_n, double b_0, arma::vec& beta, arma::vec& se);
RcppExport SEXP _bcb_lm_nig(SEXP XtySEXP, SEXP m_0SEXP, SEXP Lambda_0SEXP, SEXP Lambda_nSEXP, SEXP V_nSEXP, SEXP ytySEXP, SEXP a_nSEXP, SEXP b_0SEXP, SEXP betaSEXP, SEXP seSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type Xty(XtySEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type m_0(m_0SEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Lambda_0(Lambda_0SEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Lambda_n(Lambda_nSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type V_n(V_nSEXP);
    Rcpp::traits::input_parameter< double >::type yty(ytySEXP);
    Rcpp::traits::input_parameter< double >::type a_n(a_nSEXP);
    Rcpp::traits::input_parameter< double >::type b_0(b_0SEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type se(seSEXP);
    lm_nig(Xty, m_0, Lambda_0, Lambda_n, V_n, yty, a_n, b_0, beta, se);
    return R_NilValue;
END_RCPP
}
// lm_nig0
void lm_nig0(arma::mat& X, arma::vec& y, arma::vec& m_0, arma::mat& Lambda_0, double a_0, double b_0, double n, arma::vec& beta, arma::vec& se);
RcppExport SEXP _bcb_lm_nig0(SEXP XSEXP, SEXP ySEXP, SEXP m_0SEXP, SEXP Lambda_0SEXP, SEXP a_0SEXP, SEXP b_0SEXP, SEXP nSEXP, SEXP betaSEXP, SEXP seSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type m_0(m_0SEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Lambda_0(Lambda_0SEXP);
    Rcpp::traits::input_parameter< double >::type a_0(a_0SEXP);
    Rcpp::traits::input_parameter< double >::type b_0(b_0SEXP);
    Rcpp::traits::input_parameter< double >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type se(seSEXP);
    lm_nig0(X, y, m_0, Lambda_0, a_0, b_0, n, beta, se);
    return R_NilValue;
END_RCPP
}
// lookup
int lookup(const arma::vec parents, const arma::mat ps_i);
RcppExport SEXP _bcb_lookup(SEXP parentsSEXP, SEXP ps_iSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec >::type parents(parentsSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type ps_i(ps_iSEXP);
    rcpp_result_gen = Rcpp::wrap(lookup(parents, ps_i));
    return rcpp_result_gen;
END_RCPP
}
// lookup_score_cpp
double lookup_score_cpp(const arma::vec parents, const arma::mat ps_i, int score_col);
RcppExport SEXP _bcb_lookup_score_cpp(SEXP parentsSEXP, SEXP ps_iSEXP, SEXP score_colSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec >::type parents(parentsSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type ps_i(ps_iSEXP);
    Rcpp::traits::input_parameter< int >::type score_col(score_colSEXP);
    rcpp_result_gen = Rcpp::wrap(lookup_score_cpp(parents, ps_i, score_col));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_bcb_obj_cpp", (DL_FUNC) &_bcb_obj_cpp, 5},
    {"_bcb_con_cpp", (DL_FUNC) &_bcb_con_cpp, 4},
    {"_bcb_beta_cpp", (DL_FUNC) &_bcb_beta_cpp, 3},
    {"_bcb_lm_cpp", (DL_FUNC) &_bcb_lm_cpp, 3},
    {"_bcb_lm_nig", (DL_FUNC) &_bcb_lm_nig, 10},
    {"_bcb_lm_nig0", (DL_FUNC) &_bcb_lm_nig0, 9},
    {"_bcb_lookup", (DL_FUNC) &_bcb_lookup, 2},
    {"_bcb_lookup_score_cpp", (DL_FUNC) &_bcb_lookup_score_cpp, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_bcb(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
