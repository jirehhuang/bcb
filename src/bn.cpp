// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>



// objective function for NlcOptim::solnl()

// [[Rcpp::export]]
double obj_cpp(arma::vec par,
               arma::uvec pt_dim,
               arma::mat x_tp_,
               arma::vec b_t,
               arma::mat a_p) {

  double t, loss = 0;
  for (int i = 0; i < pt_dim(0); i++){
    t = 0;
    for (int j = 0; j < pt_dim(1); j++){
      t += a_p(i, j) * x_tp_(i, j) * par(i) * par(j + pt_dim(0));
    }
    loss += std::abs(t - b_t(i));
  }
  return loss;
}



// portion of constraint function for NlcOptim::solnl()

// [[Rcpp::export]]
Rcpp::List con_cpp(Rcpp::NumericVector par,
                   arma::uvec pt_dim,
                   arma::mat x_tp_,
                   arma::vec b_t) {

  Rcpp::NumericVector ceq(pt_dim(1));
  // arma::vec ceq(pt_dim(1), arma::fill::zeros);
  double t, loss = 0;
  for (int j = 0; j < pt_dim(1); j++){
    t = 0;
    for (int i = 0; i < pt_dim(0); i++){
      ceq(j) += x_tp_(i, j) * par(i) * par(j + pt_dim(0));
    }
    ceq(j) -= 1;
  }
  return Rcpp::List::create(Rcpp::_["ceq"] = ceq, Rcpp::_["c"] = -par);
}
