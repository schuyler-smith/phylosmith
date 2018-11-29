/*
 * 
 *
 * 
 *      Author: Schuyler Smith
 */
#include <iostream>
#include <string>
#include <regex>
#include <RcppArmadillo.h>
#include <RcppParallel.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppParallel)]]

using namespace std;

// [[Rcpp::export]]

Rcpp::DataFrame match_sequences(Rcpp::NumericMatrix short_input, Rcpp::NumericMatrix long_input){

    vector<string> shorter_seqs = Rcpp::as<vector<string>>(rownames(clone(short_input)));
    vector<string> longer_seqs = Rcpp::as<vector<string>>(rownames(clone(long_input)));
    int s_size = shorter_seqs.size();
    int l_size = longer_seqs.size();
    vector<pair<double,double>> match_seqs;
    vector<pair<double,double>> match_dups;

    #pragma omp parallel for
    for(int i=0; i<s_size; ++i){ //loop through each sequence of set with shorter read lengths
        string s_seq (shorter_seqs[i]);
        int ind = -1; //assign ind to not be within the dataframe
        for(int j=0; j<l_size; ++j){ //loop through other dataset
            string l_seq (longer_seqs[j]);
            if (l_seq != s_seq){
                if (l_seq.find(s_seq) != std::string::npos){ //if seq 1 is contained in seq 2
                    #pragma omp critical
                    if (ind == -1){
                        match_seqs.push_back(pair<double,double>(i,j));
                        ind = j;
                    } else {
                        match_dups.push_back(pair<double,double>(ind,j));
                    }
                }
            } else {
                if (ind == -1){
                    ind = j;
                } else {
                    #pragma omp critical
                    match_dups.push_back(pair<double,double>(ind,j));
                }
            }
        }
    }
// replace match of longer seqs with shorter seqs
    int match_size = match_seqs.size();
    for(int r=0; r<match_size; ++r){
        longer_seqs[match_seqs[r].second] = shorter_seqs[match_seqs[r].first];
    }
// create order of duplicates in descending index (as index changes as rows will be removed)
    int dup_size = match_dups.size();
    vector<double> duplicates;
    for(int r=0; r<dup_size; ++r){
        duplicates.push_back(match_dups[r].second);
    }
    arma::mat M(duplicates);
    arma::uvec dup_indices = sort_index(M, "descend");
    
// remove data rows and duplicate taxa names
    Rcpp::StringVector long_seqs = Rcpp::wrap(longer_seqs);
    arma::mat l_input = Rcpp::as<arma::mat>(clone(long_input));
    for(int r=0; r<dup_size; ++r){ 
        int dup_index (dup_indices[r]);
        l_input.row(match_dups[dup_index].first) = l_input.row(match_dups[dup_index].first) + l_input.row(match_dups[dup_index].second); //add values to original
        l_input.shed_row(match_dups[dup_index].second);
        long_seqs.erase(match_dups[dup_index].second);
    }
//convert back to rcpp matrix to return to R
    Rcpp::NumericMatrix x = Rcpp::wrap(l_input); 
        rownames(x) = long_seqs;
        colnames(x) = colnames(long_input);
    
    return x;
}



