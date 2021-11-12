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
