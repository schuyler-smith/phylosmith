/*
 *
 *	Author: Schuyler Smith
 *
 */
#include "significance.h"

#include <cmath>
#include <string>
#include <vector>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(RcppProgress)]]
#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include <progress.hpp>

#ifdef _OPENMP
  #include <omp.h>
#endif


//' @author Schuyler D. Smith
//' @title Assign rank values to a matrix.
//' @description Arranges the co-occurence table so that taxa of interest are
//' on the left.
//' @param count_matrix A \code{data.frame}
//' @return A rank matrix
// [[Rcpp::export]]
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

//' @author Schuyler D. Smith
//' @title Co-occurrence calculation
//' @description Calculate the pair-wise Spearman rank correlation.
//' @usage Correlation(X, Y,
//' cor_coef_cutoff, p_cutoff, method, ncores)
//' @param X An \code{otu_table} in the format from
//' \code{\link[phyloseq:otu_table]{phyloseq}}
//' @param Y An \code{otu_table} in the format from
//' \code{\link[phyloseq:otu_table]{phyloseq}}
//' @param cor_coef_cutoff \code{double} representing the minimum \code{rho-value}
//' accepted for the correlation to be returned.
//' @param p_cutoff \code{double} representing the maximum \code{p-value}
//' accepted for the correlation to be returned.
//' @param method Pearson, Spearman
//' @param ncores An \code{int} for how many cores to use to multithread the
//' calculations.
//' @return A \code{data.frame} with treatment, otu_1, otu_2, rho, p values.
//' @seealso \code{\link{co_occurrence}}
// [[Rcpp::export]]
Rcpp::DataFrame Correlation(
    Rcpp::NumericMatrix X,
    Rcpp::NumericMatrix Y = Rcpp::NumericMatrix(1),
    const double cor_coef_cutoff = 0,
    const double p_cutoff = 1,
    const std::string method = "pearson",
    const int ncores = 1
){
  std::vector<int> X_names;
  std::vector<int> Y_names;
  std::vector<double> p_values;
  std::vector<double> cor_coef_values;

  size_t N_X = X.nrow();
  size_t df = X.ncol() - 2;

  arma::mat X_comp;
  const std::string pearson = "pearson";
  const std::string spearman = "spearman";
  const std::string kendall = "kendall";
  if(method != pearson){
    X_comp = assign_rank(X);
  } else {
    X_comp = Rcpp::as<arma::mat>(clone(X));
  }
  if(Y.nrow() > 1){
    arma::mat Y_comp;
    size_t N_Y = Y.ncol();
    if(method != pearson){
      Y_comp = assign_rank(Rcpp::transpose(Y));
    } else {
      Y_comp = Rcpp::as<arma::mat>(clone(Y)).t();
    }
    Progress p(1, false);
    #ifdef _OPENMP
      #pragma omp parallel for num_threads(ncores)
    #endif
    for(size_t i=0; i<N_X; ++i){
      if(!Progress::check_abort()){
        arma::rowvec X_values = X_comp.row(i);
        for(size_t j=0; j<N_Y; ++j){
          double cor_coef;
          double p_val;
          arma::rowvec Y_values = Y_comp.row(j);
          if(arma::sum(X_values) > 0 && arma::sum(Y_values) > 0){
            std::vector<double> X_i = arma::conv_to<std::vector <double> >::from(X_values);
            std::vector<double> Y_j = arma::conv_to<std::vector <double> >::from(Y_values);
            if(method != kendall){
              cor_coef = pearsoncoeff(X_i, Y_j);
            } else {
              cor_coef = kendall_tau(X_i, Y_j);
            }
            double t = t_statistic(cor_coef, df);
            p_val = pvalue(t, df);
          } else {
            cor_coef = 0;
            p_val = 1;
          }
          if((cor_coef >= cor_coef_cutoff || cor_coef <= -cor_coef_cutoff) && p_val <= p_cutoff){ // these pushbacks takes the longest amount of time.. i think because of how memory is allocated, may need to look for more optimal method
            #ifdef _OPENMP
              #pragma omp critical
            #endif
            {
              p_values.push_back(p_val);
              cor_coef_values.push_back(cor_coef);
              X_names.push_back(i+1);
              Y_names.push_back(j+1);
            }
          }
        }
      }}
  } else {
    Progress p(1, false);
    #ifdef _OPENMP
      #pragma omp parallel for num_threads(ncores)
    #endif
    for(size_t i=0; i<N_X-1; ++i){
      if(!Progress::check_abort()){
        arma::rowvec X_values = X_comp.row(i);
        for(size_t j=i+1; j<N_X; ++j){
          double cor_coef;
          double p_val;
          arma::rowvec Y_values = X_comp.row(j);
          if(arma::sum(X_values) > 0 && arma::sum(Y_values) > 0){
            std::vector<double> X_i = arma::conv_to<std::vector <double> >::from(X_values);
            std::vector<double> Y_j = arma::conv_to<std::vector <double> >::from(Y_values);
            if(method != kendall){
              cor_coef = pearsoncoeff(X_i, Y_j);
            } else {
              cor_coef = kendall_tau(X_i, Y_j);
            }
            double t = t_statistic(cor_coef, df);
            p_val = pvalue(t, df);
          } else {
            cor_coef = 0;
            p_val = 1;
          }
          if((cor_coef >= cor_coef_cutoff || cor_coef <= -cor_coef_cutoff) && p_val <= p_cutoff){ // these pushbacks takes the longest amount of time.. i think because of how memory is allocated, may need to look for more optimal method
            #ifdef _OPENMP
              #pragma omp critical
            #endif
            {
              p_values.push_back(p_val);
              cor_coef_values.push_back(cor_coef);
              X_names.push_back(i+1);
              Y_names.push_back(j+1);
            }
          }
        }
      }}
  }
  return Rcpp::DataFrame::create(
    Rcpp::Named("X") = X_names,
    Rcpp::Named("Y") = Y_names,
    Rcpp::Named("rho") = cor_coef_values,
    Rcpp::Named("p") = p_values
  );
}

//' @author Schuyler D. Smith
//' @title Co-occurrence rho calculations
//' @description Calculates the pair-wise Spearman rank correlation without
//' testing for significance.
//' @usage permute_rho_Rcpp(count_matrix, permuted_matrix, method, ncores)
//' @param count_matrix An \code{otu_table} in the format from
//' \code{\link[phyloseq:otu_table]{phyloseq}}
//' @param permuted_matrix An \code{otu_table} in the format from
//' \code{\link[phyloseq:otu_table]{phyloseq}}
//' @param method pearson, spearman
//' @param ncores An \code{int} for how many cores to use to multithread the
//' calculations.
//' @return A \code{vector} with rho values for each pair-wise correlation.
//' @seealso \code{\link{permute_rho}}
// [[Rcpp::export]]
Rcpp::DataFrame permute_rho_Rcpp(
    Rcpp::NumericMatrix count_matrix,
    Rcpp::NumericMatrix permuted_matrix,
    const std::string method = "pearson",
    const int ncores = 1){

  size_t N = count_matrix.nrow();
  std::vector<double> cor_coef_values;

  arma::mat comparison_matrix;
  arma::mat permuted_comparison_matrix;
  const std::string pearson = "pearson";
  const std::string spearman = "spearman";
  const std::string kendall = "kendall";
  if(method != pearson){
    comparison_matrix = assign_rank(count_matrix);
    permuted_comparison_matrix = assign_rank(permuted_matrix);
  } else {
    comparison_matrix = Rcpp::as<arma::mat>(clone(count_matrix));
    permuted_comparison_matrix = Rcpp::as<arma::mat>(clone(permuted_matrix));
  }

  Progress p(1, false);
  #ifdef _OPENMP
    #pragma omp parallel for num_threads(ncores)
  #endif
  for(size_t i=0; i<N-1; ++i){
    if(!Progress::check_abort()){
      arma::rowvec X_values = comparison_matrix.row(i);
      for(size_t j=i+1; j<N; ++j){
        double cor_coef;
        arma::rowvec Y_values = permuted_comparison_matrix.row(j);
        if(arma::sum(X_values) > 0 && arma::sum(Y_values) > 0){
          std::vector<double> X = arma::conv_to<std::vector <double> >::from(X_values);
          std::vector<double> Y = arma::conv_to<std::vector <double> >::from(Y_values);
          if(method != kendall){
            cor_coef = pearsoncoeff(X, Y);
          } else {
            cor_coef = kendall_tau(X, Y);
          }
        } else {
          cor_coef = 0;
        }
        #ifdef _OPENMP
          #pragma omp critical
        #endif
        {
          cor_coef_values.push_back(cor_coef);
        }
      }}
  }
  return Rcpp::DataFrame::create(
    Rcpp::Named("rho") = cor_coef_values
  );
}

//' @author Schuyler D. Smith
//' @title Arrange co-occurence table
//' @description Arranges the co-occurence table so that taxa of interest are
//' on the left.
//' @param co_occurrence_table A \code{data.frame} in the format from
//' \code{\link{co_occurrence}}.
//' @param taxa_of_interest A \code{vector} containing names of taxa to be
//' found in either OTU_1 or OTU_2.
//' @return A \code{data.frame} with treatment, otu_1, otu_2, rho, p values.
//' @seealso \code{\link{co_occurrence}}
// [[Rcpp::export]]
Rcpp::DataFrame arrange_co_occurrence_table(Rcpp::DataFrame co_occurrence_table,
	Rcpp::CharacterVector taxa_of_interest){

    Rcpp::CharacterVector treatment = co_occurrence_table[0];
    Rcpp::CharacterVector taxa_1 = co_occurrence_table[1];
    Rcpp::CharacterVector taxa_2 = co_occurrence_table[2];
    Rcpp::NumericVector rho_values = co_occurrence_table[3];
    Rcpp::NumericVector p_values = co_occurrence_table[4];
    std::vector<std::string> headers = co_occurrence_table.names();
    size_t n_pairs = co_occurrence_table.nrow();
    size_t n_taxa = taxa_of_interest.size();

    #pragma omp parallel for
    for(int row=0; row<n_pairs; ++row){
        std::string t_1 = Rcpp::as<std::string>(taxa_1[row]);
        std::string t_2 = Rcpp::as<std::string>(taxa_2[row]);
        bool t_1_in = FALSE;
        bool t_2_in = FALSE;
        for(size_t taxa=0; taxa<n_taxa; ++taxa){
            if (t_1.find(taxa_of_interest[taxa]) != std::string::npos){
                t_1_in = TRUE;
            }
            if (t_2.find(taxa_of_interest[taxa]) != std::string::npos){
                t_2_in = TRUE;
            }
        }
        if(t_1_in == TRUE && t_2_in == FALSE){
            continue;
        } else if(t_1_in == FALSE && t_2_in == TRUE){
            taxa_1[row] = t_2;
            taxa_2[row] = t_1;
        } else if(t_1_in == TRUE && t_2_in == TRUE){
            #pragma omp critical
            {
                treatment.push_back(treatment[row]);
                taxa_1.push_back(taxa_1[row]);
                taxa_2.push_back(taxa_2[row]);
                rho_values.push_back(rho_values[row]);
                p_values.push_back(p_values[row]);
            }
        } else {
            continue;
        }
    }

	return Rcpp::DataFrame::create(
	    Rcpp::Named(headers[0]) = treatment,
	    Rcpp::Named(headers[1]) = taxa_1,
	    Rcpp::Named(headers[2]) = taxa_2,
	    Rcpp::Named(headers[3]) = rho_values,
	    Rcpp::Named(headers[4]) = p_values
	    );
}


// CHECK FOR TIED RANKS

// for(int trt=0; trt<n_treatments; ++trt){
// 	has_ties = false;
// 	arma::uvec treatment_columns = Rcpp::as<arma::uvec>(treatment_indices[trt]);
// 	arma::mat treatment_matrix = rank_table.cols(treatment_columns);
// 	map<double, int> countMap;
// 		arma::rowvec ranks = treatment_matrix.row(8);
// 		for(auto & elem : ranks){
// 			auto result = countMap.insert(pair<double, int>(elem, 1));
// 			if(result.second == false){result.first->second++;}
// 		}
// 		for(auto & elem : countMap){if(elem.first != 0 && elem.second > 1){has_ties = true;}}
// }
