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

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppParallel)]]

#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace std;

// [[Rcpp::export]]

Rcpp::DataFrame arrange_co_occurrence_table(Rcpp::DataFrame co_occurrence_table, Rcpp::CharacterVector taxa_of_interest){

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



