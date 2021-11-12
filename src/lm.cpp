// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;



// [[Rcpp::export]]
void lm_cpp(arma::mat& X,
            arma::vec& y,
            arma::vec& beta,
            arma::vec& se){

  // coefficients
  arma::vec betas = arma::solve(X, y);

  // standard errors
  arma::vec resid = y - X * betas;
  arma::vec stderrs = arma::sqrt(
    arma::dot(resid, resid) / (X.n_rows - X.n_cols - 1) *
      arma::diagvec(arma::inv(arma::trans(X) * X))
  );

  beta(0) = betas(0);
  se(0) = stderrs(0);
}



// [[Rcpp::export]]
void lm_nig(arma::vec& Xty,
            arma::vec& m_0,
            arma::mat& Lambda_0,
            arma::mat& Lambda_n,
            arma::mat& V_n,
            double yty,
            double a_n,
            double b_0,
            arma::vec& beta,
            arma::vec& se){

  // coefficients
  arma::vec m_n = V_n * (Xty + Lambda_0 * m_0);

  double b_n = b_0 + (yty +
                      arma::dot(m_0, Lambda_0 * m_0) -
                      arma::dot(m_n, Lambda_n * m_n)) / 2;

  beta(0) = m_n(0);
  se(0) = b_n / a_n * V_n(0, 0);  // standard error
}



// [[Rcpp::export]]
void lm_nig0(arma::mat& X,
             arma::vec& y,
             arma::vec& m_0,
             arma::mat& Lambda_0,
             double a_0,
             double b_0,
             double n,
             arma::vec& beta,
             arma::vec& se){

  arma::mat Lambda_n = arma::trans(X) * X + Lambda_0;
  arma::mat V_n = Lambda_n.i();

  // coefficients
  arma::vec m_n = V_n * (arma::trans(X) * y + Lambda_0 * m_0);

  double a_n = a_0 + n / 2;
  double b_n = b_0 + (arma::dot(y, y) +
                      arma::dot(m_0, Lambda_0 * m_0) -
                      arma::dot(m_n, Lambda_n * m_n)) / 2;

  beta(0) = m_n(0);
  se(0) = b_n / a_n * V_n(0, 0);  // standard error

  // return Rcpp::List::create(Named("m_n") = m_n,
  //                           Named("Lambda_n") = Lambda_n,
  //                           Named("a_n") = a_n,
  //                           Named("b_n") = b_n,
  //                           Named("V_n") = V_n,
  //                           Named("Xty") = arma::trans(X) * y,
  //                           Named("Lm_0") = Lambda_0 * m_0,
  //                           Named("yTy") = arma::dot(y, y),
  //                           Named("mLm_0") = arma::dot(m_0, Lambda_0 * m_0),
  //                           Named("mLm_n") = arma::dot(m_n, Lambda_n * m_n));
}
