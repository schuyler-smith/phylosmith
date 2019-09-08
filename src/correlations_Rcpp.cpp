/*
 *
 *	Author: Schuyler Smith
 *
 */
#include <cmath>
#include <vector>
#include <string>
#include <regex>
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(RcppProgress)]]
#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include <progress.hpp>

#ifdef _OPENMP
  #include <omp.h>
#endif

#include "significance.h"
#include "rank.h"

using namespace std;

//' @author Schuyler D. Smith
//' @title Co-occurrence calculation
//' @description Calculate the pair-wise Spearman rank correlation.
//' @usage Correlation(otu_table, treatment_indices,
//' rho_cutoff, p_cutoff, method, ncores)
//' @param otu_table An \code{otu_table} in the format from
//' \code{\link[phyloseq:otu_table]{phyloseq}}
//' @param treatment_indices A \code{list} with c++ indices for the
//' \code{treatment_names} corresponding to which treatment each column in the
//' \code{otu_table} belongs to.
//' @param treatment_names A \code{Vector} containing the treatment names
//' corresponding to the \code{treatment_indices}.
//' @param rho_cutoff \code{double} representing the minimum \code{rho-value}
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
    Rcpp::NumericMatrix count_matrix,
    const double cor_coef_cutoff = 0,
    const double p_cutoff = 1,
    const std::string method = "pearson",
    const int ncores = 1
){
  vector<string> X_names;
  vector<string> Y_names;
  vector<double> p_values;
  vector<double> cor_coef_values;

  vector<string> names = Rcpp::as<vector <string> >(rownames(count_matrix));
  size_t N = count_matrix.nrow();

  arma::mat comparison_matrix;
  const string pearson = "pearson";
  const string spearman = "spearman";
  const string kendall = "kendall";
  if(method != pearson){
    comparison_matrix = assign_rank(count_matrix);
  } else {
    comparison_matrix = Rcpp::as<arma::mat>(clone(count_matrix));
  }

  Progress p(1, false);
  int n_samples = comparison_matrix.n_cols;
  #ifdef _OPENMP
    #pragma omp parallel for num_threads(ncores)
  #endif
  for(size_t subject_1=0; subject_1<N-1; ++subject_1){
    if(!Progress::check_abort()){
      arma::rowvec X_values = comparison_matrix.row(subject_1);
      for(size_t subject_2=subject_1+1; subject_2<N; ++subject_2){
        double cor_coef;
        double p_val;
        arma::rowvec Y_values = comparison_matrix.row(subject_2);
        if(arma::sum(X_values) > 0 && arma::sum(Y_values) > 0){
          vector<double> X = arma::conv_to<vector <double> >::from(X_values);
          vector<double> Y = arma::conv_to<vector <double> >::from(Y_values);
          if(method != kendall){
            cor_coef = pearsoncoeff(X, Y);
          } else {
            cor_coef = kendall_tau(X, Y);
          }
          double df = n_samples - 2;
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
            X_names.push_back(names[subject_1]);
            Y_names.push_back(names[subject_2]);
          }
        }
      }
    }}

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
//' @usage permute_rho_Rcpp(otu_table, treatment_indices,
//' treatment_names, ncores)
//' @param otu_table An \code{otu_table} in the format from
//' \code{\link[phyloseq:otu_table]{phyloseq}}
//' @param treatment_indices A \code{list} with c++ indices for the
//' \code{treatment_names} corresponding to which treatment each column in the
//' \code{otu_table} belongs to.
//' @param treatment_names A \code{Vector} containing the treatment names
//' corresponding to the \code{treatment_indices}.
//' @param method Pearson, Spearman
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
  vector<double> cor_coef_values;

  arma::mat comparison_matrix;
  arma::mat permuted_comparison_matrix;
  const string pearson = "pearson";
  const string spearman = "spearman";
  const string kendall = "kendall";
  if(method != pearson){
    comparison_matrix = assign_rank(count_matrix);
    permuted_comparison_matrix = assign_rank(permuted_matrix);
  } else {
    comparison_matrix = Rcpp::as<arma::mat>(clone(count_matrix));
    permuted_comparison_matrix = Rcpp::as<arma::mat>(clone(permuted_matrix));
  }

  Progress p(1, false);
  int n_samples = comparison_matrix.n_cols;
  #ifdef _OPENMP
    #pragma omp parallel for num_threads(ncores)
  #endif
  for(size_t subject_1=0; subject_1<N-1; ++subject_1){
    if(!Progress::check_abort()){
      arma::rowvec X_values = comparison_matrix.row(subject_1);
      for(size_t subject_2=subject_1+1; subject_2<N; ++subject_2){
        double cor_coef;
        arma::rowvec Y_values = permuted_comparison_matrix.row(subject_2);
        if(arma::sum(X_values) > 0 && arma::sum(Y_values) > 0){
          vector<double> X = arma::conv_to<vector <double> >::from(X_values);
          vector<double> Y = arma::conv_to<vector <double> >::from(Y_values);
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
    vector<string> headers = co_occurrence_table.names();
    int n_pairs = co_occurrence_table.nrow();
    int n_taxa = taxa_of_interest.size();

    #pragma omp parallel for
    for(int row=0; row<n_pairs; ++row){
        string t_1 = Rcpp::as<string>(taxa_1[row]);
        string t_2 = Rcpp::as<string>(taxa_2[row]);
        bool t_1_in = FALSE;
        bool t_2_in = FALSE;
        for(int taxa=0; taxa<n_taxa; ++taxa){
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
