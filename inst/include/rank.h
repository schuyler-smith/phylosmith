#ifndef RANK_H
#define RANK_H

#include "significance.h"

#include <cmath>
#include <vector>
#include <string>
#include <regex>
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(RcppProgress)]]
#include <RcppArmadillo.h>

using namespace std;

arma::mat assign_rank(Rcpp::NumericMatrix count_matrix
){
  arma::mat rank_table = Rcpp::as<arma::mat>(clone(count_matrix));
  size_t n_X = count_matrix.nrow();

  arma::uvec ordered_indices;
  int rank;
  int ties = 1;

  size_t n_samples = rank_table.n_cols;
  for(size_t X=0; X<n_X; ++X){
    arma::rowvec rank_vector = arma::zeros<arma::rowvec>(n_samples);
    arma::rowvec X_abund = rank_table.row(X);
    if(arma::sum(X_abund) > 0){
      ordered_indices = sort_index(X_abund, "descend"); // sort indices from largest to smallest
      for(size_t rank_index=0; rank_index<n_samples; ++rank_index){ // indices are confusing here, likely source of any problem
        int X_index = ordered_indices[rank_index];
        rank = rank_index+1;
        if(rank_index < (n_samples-1) && X_abund[X_index] == X_abund[ordered_indices[rank]]){ // especially here
          ++ties;
          rank_vector[X_index] = rank;
          continue; // if the next sample has an equal value, go to next iteration
        }
        rank_vector[X_index] = rank;
        if(ties > 1){
          arma::uvec tied_indices = arma::find(rank_vector > rank-ties && rank_vector <= rank);
          rank_vector.elem(tied_indices).fill(arma::sum(rank_vector.elem(tied_indices)) / ties);
          ties = 1; // sum all equal values' ranks and assign them as the mean
        }
      }
    }
    rank_table.row(X) = rank_vector; // put ranks into matrix
  }
  return rank_table;
}

double kendall_tau(vector<double> X,
                   vector<double> Y){
  size_t N = X.size();
  vector<double> sorted_Y;
  sorted_Y.resize(N);
  for (size_t i = 0; i<N; i++){
    sorted_Y[X[i]-1] = Y[i];
  }
  vector<double> concordant;
  vector<double> discordant;
  for (size_t i = 0; i<N-1; i++){
    for (size_t j = i+1; j<N; j++){
      discordant.push_back(sorted_Y[i] > sorted_Y[j]);
      concordant.push_back(sorted_Y[i] < sorted_Y[j]);
    }
  }
  return (sum(concordant) - sum(discordant))/(sum(concordant) + sum(discordant));
}



#endif /* RANK_H */
