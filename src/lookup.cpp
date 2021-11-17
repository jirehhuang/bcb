// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;



// [[Rcpp::export]]
double lookup(const arma::vec parents,
              const arma::mat cache){

  int n_parents = parents.size();

  // no parents
  if (n_parents == 0)
    return 1;

  int j = 0; // iterate over elements of parents

  // iterate over each parent set
  int i = 0;
  while (i < cache.n_rows){

    // skip if last parent element is NA
    if (!arma::is_finite(cache(i, n_parents - 1))){

      i++;
      continue;
    }

    if (cache(i, j) == parents(j)){

      // found exact row
      if (j == (n_parents - 1)) break;

      j++;  // search next element
      continue;
    }
    i++;
  }

  // if successful, return score
  if (i < cache.n_rows)
    return i + 1;

  return 0;
}



// [[Rcpp::export]]
double lookup_score_cpp(const arma::vec parents,
                        const arma::mat cache){

  int i = lookup(parents, cache) - 1L;

  return cache(i, cache.n_cols - 1);

  // int n_parents = parents.size();
  //
  // // no parents
  // if (n_parents == 0)
  //   return cache(0, cache.n_cols - 1);
  //
  // int j = 0; // iterate over elements of parents
  //
  // // iterate over each parent set
  // int i = 0;
  // while (i < cache.n_rows){
  //
  //   // skip if last parent element is NA
  //   if (!arma::is_finite(cache(i, n_parents - 1))){
  //
  //     i++;
  //     continue;
  //   }
  //
  //   if (cache(i, j) == parents(j)){
  //
  //     // found exact row
  //     if (j == (n_parents - 1)) break;
  //
  //     j++;  // search next element
  //     continue;
  //   }
  //   i++;
  // }
  //
  // // if successful, return score
  // if (i < cache.n_rows)
  //   return cache(i, cache.n_cols - 1);
  //
  // return 0;
}
