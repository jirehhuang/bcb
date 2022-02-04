// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;



// [[Rcpp::export]]
void beta_cpp(arma::mat& X,
              arma::vec& y,
              arma::vec& beta){

  // coefficients
  arma::vec betas = arma::solve(X, y);
  for (int i = 0; i < beta.size(); i++){

    beta(i) = betas(i);
  }
}



// [[Rcpp::export]]
void lm_cpp(arma::mat& X,
            arma::vec& y,
            arma::vec& values){

  // coefficients
  arma::vec betas = arma::solve(X, y);
  values(0) = betas(0);

  // standard errors
  arma::vec resid = y - X * betas;
  arma::mat inv_XtX = arma::inv(X.t() * X);
  values(2) = arma::dot(resid, resid);  // sum of squared residuals
  values(3) = 1 / inv_XtX(0, 0);  // XtX
  values(1) = std::sqrt(values(2) / (X.n_rows - X.n_cols) / values(3));

  // solve for SE_bda = SE_nig
  resid = y - X.col(0) * betas(0);
  resid = resid - arma::mean(resid);
  double values_2 = arma::dot(resid, resid);  // update to estimate interventional sd
  values(4) = (values_2 / (X.n_rows - X.n_cols - 1)) *  // 2b_0 / (2a_0)
    (1 / values(2)) *  // 1 / RSS
    (X.n_rows - X.n_cols) * values(3);  // (n - k) * (1 / inv_XtX)
  values(2) = values_2;
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

  // for (int i = 0; i < beta.size(); i++){
  //
  //   beta(i) = m_n(i);
  //   se(i) = std::sqrt(b_n / a_n * V_n(i, i));  // standard error
  // }
  beta(0) = m_n(0);
  se(0) = std::sqrt(b_n / a_n * V_n(0));
  // beta(beta.size() - 1) = m_n(m_n.size() - 1);
  // se(se.size() - 1) = std::sqrt(b_n / a_n * V_n(V_n.size() - 1));
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
