/*
 * 
 *	Author: Schuyler Smith
 *  significance.h written by Jin Choi
 *      
 */
#include "significance.h"
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(RcppProgress)]]
#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include <progress.hpp>

#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace std;

// [[Rcpp::export]]

std::vector<double> FastCoOccur_rho_Rcpp(
	Rcpp::NumericMatrix otu_table, 
	Rcpp::List treatment_indices, 
	Rcpp::StringVector treatment_names,
	const int ncores){

	arma::mat rank_table = Rcpp::as<arma::mat>(clone(otu_table));
	int n_treatments = treatment_names.length();
	int n_taxa = otu_table.nrow();

	arma::uvec ordered_indices;
	int rank;
	int ties = 1;

	// return vectors
	vector<double> rho_values;

	// bool has_ties;
	// creating ranks
	for(int trt=0; trt<n_treatments; ++trt){ // loop through each treatment
		arma::uvec treatment_columns = Rcpp::as<arma::uvec>(treatment_indices[trt]); // vector of elements for each sample in this treatment
		arma::mat treatment_matrix = rank_table.cols(treatment_columns); // subset the matrix to just those samples
		int n_samples = treatment_columns.size();
		for(int taxa=0; taxa<n_taxa; ++taxa){ // loop through all the taxa
			arma::rowvec rank_vector = arma::zeros<arma::rowvec>(treatment_columns.size());
			arma::rowvec taxa_abund = treatment_matrix.row(taxa);
			if(arma::sum(taxa_abund) > 0){
				ordered_indices = sort_index(taxa_abund, "descend"); // sort indices from largest to smallest
				for(int rank_index=0; rank_index<n_samples; ++rank_index){ // indices are confusing here, likely source of any problem
					int taxa_index = ordered_indices[rank_index];
					rank = rank_index+1;
					if(rank_index < (n_samples-1) && taxa_abund[taxa_index] == taxa_abund[ordered_indices[rank]]){ // especially here
						++ties;
						rank_vector[taxa_index] = rank;
						continue; // if the next sample has an equal value, go to next iteration
					}
					rank_vector[taxa_index] = rank;
					if(ties > 1){
						arma::uvec tied_indices = arma::find(rank_vector > rank-ties && rank_vector <= rank);
						rank_vector.elem(tied_indices).fill(arma::sum(rank_vector.elem(tied_indices)) / ties);
						ties = 1; // sum all equal values' ranks and assign them as the mean
					}
				}
			}
			treatment_matrix.row(taxa) = rank_vector; // put ranks into matrix
		}
		rank_table.cols(treatment_columns) = treatment_matrix; // put treatment back into whole table
	}

	Progress p(n_treatments, false);
	for(int trt=0; trt<n_treatments; ++trt){
		// has_ties = false;
		arma::uvec treatment_columns = Rcpp::as<arma::uvec>(treatment_indices[trt]);
		arma::mat treatment_matrix = rank_table.cols(treatment_columns);
		int n_samples = treatment_columns.size();
		#ifdef _OPENMP
			#pragma omp parallel for num_threads(ncores)
		#endif
		for(int taxa1=0; taxa1<n_taxa-1; ++taxa1){
			if(!Progress::check_abort()){
			arma::rowvec taxa1_ranks = treatment_matrix.row(taxa1);
			for(int taxa2=taxa1+1; taxa2<n_taxa; ++taxa2){
				double rho;
				arma::rowvec taxa2_ranks = treatment_matrix.row(taxa2);
				if(arma::sum(taxa1_ranks) > 0 && arma::sum(taxa2_ranks) > 0){
					vector<double> X = arma::conv_to<vector <double> >::from(taxa1_ranks);
					vector<double> Y = arma::conv_to<vector <double> >::from(taxa2_ranks);
					rho = pearsoncoeff(X, Y);
				} else {
					rho = 0; 
				}
				#ifdef _OPENMP
					#pragma omp critical
				#endif
				{
					rho_values.push_back(rho);
				}
			}	
		}}
	}

return(rho_values);
}

