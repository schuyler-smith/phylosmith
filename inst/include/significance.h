/*
 * 
 *  Author: Jin Choi
 *  Maintainer: Schuyler Smith
 *      
 */

#ifndef SIGNIFICANCE_H
#define SIGNIFICANCE_H

#include <vector>

using namespace std;

double sum(vector<double> a);
double mean(vector<double> a);
double sqsum(vector<double> a);
double stdev(vector<double> nums);
vector<double> operator-(vector<double> a, double b);
vector<double> operator*(vector<double> a, vector<double> b);
double pearsoncoeff(vector<double> X, vector<double> Y);
double betai(double a, double b, double x);
double betacf(double a, double b, double x);
double gammln(double xx);
double pvalue( double t, double df );

#endif /* SIGNIFICANCE_H */