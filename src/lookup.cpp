// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;



// [[Rcpp::export]]
int lookup(const arma::vec parents,
           const arma::mat ps_i){

  int n_parents = parents.size();

  // no parents
  if (n_parents == 0)
    return 1;

  int j = 0; // iterate over elements of parents

  // iterate over each parent set
  int i = 0;
  while (i < ps_i.n_rows){

    // skip if last parent element is NA
    if (!arma::is_finite(ps_i(i, n_parents - 1))){

      i++;
      continue;
    }

    if (ps_i(i, j) == parents(j)){

      // found exact row
      if (j == (n_parents - 1)) break;

      j++;  // search next element
      continue;
    }
    i++;
  }

  // if successful, return score
  if (i < ps_i.n_rows)
    return i + 1;

  return 0;
}



// [[Rcpp::export]]
double lookup_score_cpp(const arma::vec parents,
                        const arma::mat ps_i,
                        int score_col = -1){  // cpp index

  int i = lookup(parents, ps_i) - 1L;

  // [default] structure of ps: parents | score | prob | ordering (last col)
  if (score_col < 0)
    score_col = ps_i.n_cols - 3;

  // structure of scores: parents | score (last col)

  return ps_i(i, score_col);
}
