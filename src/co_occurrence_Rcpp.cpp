/*
 * 
 *	Author: Schuyler Smith
 *  significance.h written by Jin Choi
 *      
 */
#include "significance.h"
#include <iostream>
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

using namespace std;

//' @author Schuyler D. Smith
//' @title Co-occurrence calculation
//' @description Calculate the pair-wise Spearman rank correlation.
//' @param otu_table An \code{otu_table} in the format from \code{\link[phyloseq:otu_table]{phyloseq}}
//' @param treatment_indices A \code{list} with c++ indices for the \code{treatment_names} corresponding to which treatment each column in the \code{otu_table} belongs to.
//' @param treatment_names A \code{Vector} containing the treatment names corresponding to the \code{treatment_indices}.
//' @param p_cutoff \code{double} representing the maximum \code{p-value} accepted for the correlation to be returned.
//' @param ncores \code{int} for how many cores to use to multithread the calculations.
//' @return A \code{data.frame} with treatment, otu_1, otu_2, rho, p values.
//' @seealso \code{\link{co_occurrence}}
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

//' @author Schuyler D. Smith
//' @title Co-occurrence rho calculations
//' @description Calculates the pair-wise Spearman rank correlation without testing for significance.
//' @param otu_table An \code{otu_table} in the format from \code{\link[phyloseq:otu_table]{phyloseq}}
//' @param treatment_indices A \code{list} with c++ indices for the \code{treatment_names} corresponding to which treatment each column in the \code{otu_table} belongs to.
//' @param treatment_names A \code{Vector} containing the treatment names corresponding to the \code{treatment_indices}.
//' @param ncores \code{int} for how many cores to use to multithread the calculations.
//' @return A \code{vector} with rho values for each pair-wise correlation.
//' @seealso \code{\link{bootstrap_rho}}
// [[Rcpp::export]]
std::vector<double> co_occurrence_rho_Rcpp(
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

//' @author Schuyler D. Smith
//' @title Arrange co-occurence table
//' @description Arranges the co-occurence table so that taxa of interest are on the left.
//' @param co_occurrence_table A \code{data.frame} in the format from \code{\link{co_occurrence}}.
//' @param taxa_of_interest A \code{vector} containing names of taxa to be found in either OTU_1 or OTU_2.
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