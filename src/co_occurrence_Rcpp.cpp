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

Rcpp::DataFrame co_occurrence_Rcpp(
	Rcpp::NumericMatrix otu_table, 
	Rcpp::List treatment_indices, 
	Rcpp::StringVector treatment_names, 
	double p_cutoff, 
	const int ncores){

	vector<string> taxa_names = Rcpp::as<vector <string> >(rownames(otu_table));
	arma::mat rank_table = Rcpp::as<arma::mat>(clone(otu_table));
	int n_treatments = treatment_names.length();
	int n_taxa = otu_table.nrow();

	arma::uvec ordered_indices;
	int rank;
	int ties = 1;

	// return vectors
	vector<string> treatments;
	vector<double> p_values;
	vector<double> rho_values;
	vector<string> taxa_1;
	vector<string> taxa_2;

	// bool has_ties;
	// creating ranks
	for(int trt=0; trt<n_treatments; ++trt){ // loop through each treatment
		arma::uvec treatment_columns = Rcpp::as<arma::uvec>(treatment_indices[trt]); // vector of elements for each sample in this treatment
		arma::mat treatment_matrix = rank_table.cols(treatment_columns); // subset the matrix to just those samples
		int n_samples = treatment_columns.size();
		// #pragma omp parallel for
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
	  		// #pragma omp parallel for
		#endif
		for(int taxa1=0; taxa1<n_taxa-1; ++taxa1){
			if(!Progress::check_abort()){
			arma::rowvec taxa1_ranks = treatment_matrix.row(taxa1);
			for(int taxa2=taxa1+1; taxa2<n_taxa; ++taxa2){
				double rho;
				double p_val;
				arma::rowvec taxa2_ranks = treatment_matrix.row(taxa2);
				if(arma::sum(taxa1_ranks) > 0 && arma::sum(taxa2_ranks) > 0){
					// may add check for duplicate ranks, if none exist, can use below formula
					// double rho = 1 - (6*arma::sum(arma::square(taxa1_ranks - taxa2_ranks))) / (n_samples*(pow(n_samples,2)-1));
					// arma::mat matrho = arma::cov(taxa1_ranks, taxa2_ranks) / (arma::stddev(taxa1_ranks)*arma::stddev(taxa1_ranks));
					// rho = arma::conv_to<double>::from(matrho);
					vector<double> X = arma::conv_to<vector <double> >::from(taxa1_ranks);
					vector<double> Y = arma::conv_to<vector <double> >::from(taxa2_ranks);
					rho = pearsoncoeff(X, Y);
					double t = rho * sqrt((n_samples-2)/(1 - rho*rho));
	  				double df = n_samples - 2;
	  				p_val = pvalue( t, df );
				} else {
					rho = 0; 
					p_val = 1;
				}
				if(p_val <= p_cutoff){ // these pushbacks takes the longest amount of time.. i think because of how memory is allocated, may need to look for more optimal method
					#ifdef _OPENMP
				  		#pragma omp critical
					#endif
					{
						treatments.push_back(Rcpp::as<string> (treatment_names[trt]));
						p_values.push_back(p_val);
						rho_values.push_back(rho);
						taxa_1.push_back(taxa_names[taxa1]);
						taxa_2.push_back(taxa_names[taxa2]);
					}
				}
			}	
		}}
	}

return Rcpp::DataFrame::create(
	Rcpp::Named("Treatment") = treatments,
	Rcpp::Named("OTU_1") = taxa_1,
	Rcpp::Named("OTU_2") = taxa_2,
	Rcpp::Named("rho") = rho_values,
	Rcpp::Named("p") = p_values
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