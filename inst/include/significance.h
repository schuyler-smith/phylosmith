#ifndef SIGNIFICANCE_H_SDS
#define SIGNIFICANCE_H_SDS

#include <vector>

double sum(std::vector<double>);

double mean(std::vector<double>);

int triangle_number(int);

std::vector<double> operator-(std::vector<double>, double);

std::vector<double> operator*(std::vector<double>, std::vector<double>);

double stdev(std::vector<double>);

double pearsoncoeff(std::vector<double>, std::vector<double>);

double kendall_tau(std::vector<double>,std::vector<double>);

double t_statistic(double, int);

double betai(double, double, double);

double betacf(double, double, double);

double gammln(double);

double pvalue(double, double);

#endif /* SIGNIFICANCE_H_SDS */
