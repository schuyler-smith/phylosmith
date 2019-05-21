/*
 *
 *	Author: Schuyler Smith
 *  significance.h written by Jin Choi
 *
 */
#include "significance.h"
#include <cmath>
#include <vector>
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
//' @usage co_occurrence_Rcpp(otu_table, treatment_indices,
//' treatment_names, rho_cutoff, p_cutoff, ncores)
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
//' @param ncores An \code{int} for how many cores to use to multithread the
//' calculations.
//' @return A \code{data.frame} with treatment, otu_1, otu_2, rho, p values.
//' @seealso \code{\link{co_occurrence}}
// [[Rcpp::export]]
Rcpp::DataFrame co_occurrence_Rcpp(
	Rcpp::NumericMatrix otu_table,
	Rcpp::List treatment_indices,
	Rcpp::StringVector treatment_names,
	double rho_cutoff,
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
				if(rho >= rho_cutoff || rho <= -rho_cutoff && p_val <= p_cutoff){ // these pushbacks takes the longest amount of time.. i think because of how memory is allocated, may need to look for more optimal method
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

//' @author Schuyler D. Smith
//' @title Co-occurrence rho calculations
//' @description Calculates the pair-wise Spearman rank correlation without
//' testing for significance.
//' @usage co_occurrence_rho_Rcpp(otu_table, treatment_indices,
//' treatment_names, ncores)
//' @param otu_table An \code{otu_table} in the format from
//' \code{\link[phyloseq:otu_table]{phyloseq}}
//' @param treatment_indices A \code{list} with c++ indices for the
//' \code{treatment_names} corresponding to which treatment each column in the
//' \code{otu_table} belongs to.
//' @param treatment_names A \code{Vector} containing the treatment names
//' corresponding to the \code{treatment_indices}.
//' @param ncores An \code{int} for how many cores to use to multithread the
//' calculations.
//' @return A \code{vector} with rho values for each pair-wise correlation.
//' @seealso \code{\link{permute_rho}}
// [[Rcpp::export]]
Rcpp::DataFrame co_occurrence_rho_Rcpp(
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
	vector<string> treatments;
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
	  		// #pragma omp parallel for
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
				treatments.push_back(Rcpp::as<string> (treatment_names[trt]));
				rho_values.push_back(rho);
				}
			}
		}}
	}
	return Rcpp::DataFrame::create(
		Rcpp::Named("Treatment") = treatments,
		Rcpp::Named("rho") = rho_values
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

double sum(vector<double> a)
{
  double s = 0;
  for (size_t i = 0; i < a.size(); i++)
  {
    s += a[i];
  }
  return s;
}

double mean(vector<double> a)
{
  return sum(a) / a.size();
}


double sqsum(vector<double> a)
{
  double s = 0;
  for (size_t i = 0; i < a.size(); i++)
  {
    s += pow(a[i], 2);
  }
  return s;
}

double stdev(vector<double> nums)
{
  double N = nums.size();
  return pow(sqsum(nums) / N - pow(sum(nums) / N, 2), 0.5);
}

vector<double> operator-(vector<double> a, double b)
{
  vector<double> retvect;
  for (size_t i = 0; i < a.size(); i++)
  {
    retvect.push_back(a[i] - b);
  }
  return retvect;
}

vector<double> operator*(vector<double> a, vector<double> b)
{
  vector<double> retvect;
  for (size_t i = 0; i < a.size() ; i++)
  {
     retvect.push_back(a[i] * b[i]);
  }
  return retvect;
}

double pearsoncoeff(vector<double> X, vector<double> Y)
//calculate pearson coefficient
{
  return sum((X - mean(X))*(Y - mean(Y))) / (X.size()*stdev(X)* stdev(Y));
}

double betai(double a, double b, double x)
// Returns the incomplete beta function Ix(a, b).
{
  double betacf(double a, double b, double x);
  double gammln(double xx);
  double bt;
  if (x == 0.0 || x == 1.0)
    bt=0.0;
  else // Factors in front of the continued fraction.
    bt=exp(gammln(a+b)-gammln(a)-gammln(b)+a*log(x)+b*log(1.0-x));
  if (x < (a+1.0)/(a+b+2.0)) // Use continued fraction directly.
    return bt*betacf(a,b,x)/a;
  else // Use continued fraction after making the sym
    return 1.0-bt*betacf(b,a,1.0-x)/b; // metry transformation.
}


#define MAXIT 1000
#define EPS 3.0e-7
#define FPMIN 1.0e-30
double betacf(double a, double b, double x)
// Used by betai: Evaluates continued fraction for incomplete beta function by modiﬁed Lentz’s method (§5.2).
{
  int m,m2;
  double aa,c,d,del,h,qab,qam,qap;
  qab=a+b; // These q’s will be used in factors that occur
  qap=a+1.0; // in the coeﬃcients (6.4.6).
  qam=a-1.0;
  c=1.0; // First step of Lentz’s method.
  d=1.0-qab*x/qap;
  if (fabs(d) < FPMIN)
    d=FPMIN;
  d=1.0/d;
  h=d;
  for (m=1;m<=MAXIT;m++) {
    m2=2*m;
    aa=m*(b-m)*x/((qam+m2)*(a+m2));
    d=1.0+aa*d; // One step (the even one) of the recurrence.
    if (fabs(d) < FPMIN) d=FPMIN;
    c=1.0+aa/c;
    if (fabs(c) < FPMIN) c=FPMIN;
    d=1.0/d;
    h *= d*c;
    aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
    d=1.0+aa*d; // Next step of the recurrence (the odd one).
    if (fabs(d) < FPMIN) d=FPMIN;
    c=1.0+aa/c;
    if (fabs(c) < FPMIN) c=FPMIN;
    d=1.0/d;
    del=d*c;
    h *= del;
    if (fabs(del-1.0) < EPS) break; // Are we done?
  }
  return h;
}

double gammln(double xx)
// Returns the value ln[Γ(xx)] for xx > 0.
{
  // Internal arithmetic will be done in double precision, a nicety that you can omit if ﬁve-ﬁgure
  // accuracy is good enough.
  double x,y,tmp,ser;
  static double cof[6]={76.18009172947146,-86.50532032941677,24.01409824083091,-1.231739572450155,0.1208650973866179e-2,-0.5395239384953e-5};
  int j;
  y=x=xx;
  tmp=x+5.5;
  tmp -= (x+0.5)*log(tmp);
  ser=1.000000000190015;
  for (j=0;j<=5;j++) ser += cof[j]/++y;
  return -tmp+log(2.5066282746310005*ser/x);
}

double pvalue( double t, double df )
// Compute p-value of t-statistic
{
  return betai(0.5*df,0.5,df/(df+t*t));
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
