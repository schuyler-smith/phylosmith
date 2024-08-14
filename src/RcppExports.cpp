// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// assign_rank
arma::mat assign_rank(Rcpp::NumericMatrix count_matrix);
RcppExport SEXP _phylosmith_assign_rank(SEXP count_matrixSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type count_matrix(count_matrixSEXP);
    rcpp_result_gen = Rcpp::wrap(assign_rank(count_matrix));
    return rcpp_result_gen;
END_RCPP
}
// Correlation
Rcpp::DataFrame Correlation(Rcpp::NumericMatrix X, Rcpp::NumericMatrix Y, const double cor_coef_cutoff, const double p_cutoff, const std::string method, const int ncores);
RcppExport SEXP _phylosmith_Correlation(SEXP XSEXP, SEXP YSEXP, SEXP cor_coef_cutoffSEXP, SEXP p_cutoffSEXP, SEXP methodSEXP, SEXP ncoresSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const double >::type cor_coef_cutoff(cor_coef_cutoffSEXP);
    Rcpp::traits::input_parameter< const double >::type p_cutoff(p_cutoffSEXP);
    Rcpp::traits::input_parameter< const std::string >::type method(methodSEXP);
    Rcpp::traits::input_parameter< const int >::type ncores(ncoresSEXP);
    rcpp_result_gen = Rcpp::wrap(Correlation(X, Y, cor_coef_cutoff, p_cutoff, method, ncores));
    return rcpp_result_gen;
END_RCPP
}
// permute_rho_Rcpp
Rcpp::DataFrame permute_rho_Rcpp(Rcpp::NumericMatrix count_matrix, Rcpp::NumericMatrix permuted_matrix, const std::string method, const int ncores);
RcppExport SEXP _phylosmith_permute_rho_Rcpp(SEXP count_matrixSEXP, SEXP permuted_matrixSEXP, SEXP methodSEXP, SEXP ncoresSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type count_matrix(count_matrixSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type permuted_matrix(permuted_matrixSEXP);
    Rcpp::traits::input_parameter< const std::string >::type method(methodSEXP);
    Rcpp::traits::input_parameter< const int >::type ncores(ncoresSEXP);
    rcpp_result_gen = Rcpp::wrap(permute_rho_Rcpp(count_matrix, permuted_matrix, method, ncores));
    return rcpp_result_gen;
END_RCPP
}
// arrange_co_occurrence_table
Rcpp::DataFrame arrange_co_occurrence_table(Rcpp::DataFrame co_occurrence_table, Rcpp::CharacterVector taxa_of_interest);
RcppExport SEXP _phylosmith_arrange_co_occurrence_table(SEXP co_occurrence_tableSEXP, SEXP taxa_of_interestSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type co_occurrence_table(co_occurrence_tableSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type taxa_of_interest(taxa_of_interestSEXP);
    rcpp_result_gen = Rcpp::wrap(arrange_co_occurrence_table(co_occurrence_table, taxa_of_interest));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_phylosmith_assign_rank", (DL_FUNC) &_phylosmith_assign_rank, 1},
    {"_phylosmith_Correlation", (DL_FUNC) &_phylosmith_Correlation, 6},
    {"_phylosmith_permute_rho_Rcpp", (DL_FUNC) &_phylosmith_permute_rho_Rcpp, 4},
    {"_phylosmith_arrange_co_occurrence_table", (DL_FUNC) &_phylosmith_arrange_co_occurrence_table, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_phylosmith(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
