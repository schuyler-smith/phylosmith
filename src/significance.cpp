#include "significance.h"

#include <cmath>
#include <vector>

double sum(std::vector<double> a){
  double s = 0;
  for (size_t i = 0; i < a.size(); i++){
    s += a[i];
  }
  return s;
}

double mean(std::vector<double> a){
  return sum(a) / a.size();
}

int triangle_number(int x){
  return (x*(x+1)/2);
}

std::vector<double> operator-(std::vector<double> a, double b){
  std::vector<double> product_vector;
  size_t N = a.size();
  product_vector.resize(N);
  for (size_t i = 0; i < N; i++){
    product_vector[i] = (a[i] - b);
  }
  return product_vector;
}

std::vector<double> operator*(std::vector<double> a, std::vector<double> b){
  std::vector<double> product_vector;
  size_t N = a.size();
  product_vector.resize(N);
  for (size_t i = 0; i<N; i++){
    product_vector[i] = (a[i] * b[i]);
  }
  return product_vector;
}

double stdev(std::vector<double> X){
  size_t N = X.size();
  double x_mean = mean(X);
  for(std::size_t i=0; i<N; ++i){
    X[i] = std::pow(X[i] - x_mean, 2);
  }
  return std::pow(sum(X)/(N), 0.5);
}

double pearsoncoeff(std::vector<double> X, std::vector<double> Y){
  return sum((X - mean(X))*(Y - mean(Y))) / (X.size()*stdev(X)* stdev(Y));
}

double kendall_tau(std::vector<double> X,
                   std::vector<double> Y){
  size_t N = X.size();
  std::vector<double> sorted_Y;
  sorted_Y.resize(N);
  for (size_t i = 0; i<N; i++){
    sorted_Y[X[i]-1] = Y[i];
  }
  std::vector<double> concordant;
  std::vector<double> discordant;
  for (size_t i = 0; i<N-1; i++){
    for (size_t j = i+1; j<N; j++){
      discordant.push_back(sorted_Y[i] > sorted_Y[j]);
      concordant.push_back(sorted_Y[i] < sorted_Y[j]);
    }
  }
  return (sum(concordant) - sum(discordant))/(sum(concordant) + sum(discordant));
}

double t_statistic(double r, int df){
  double t = sqrt(df) * r / sqrt(1 - (r*r));
  return(t);
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



double betacf(double a, double b, double x)
  // Used by betai: Evaluates continued fraction for incomplete beta function by modiﬁed Lentz’s method (§5.2).
{
  constexpr double MAXIT { 1000 };
  constexpr double EPS { 3.0e-7 };
  constexpr double FPMIN { 1.0e-30 };
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
  double x, y, tmp,ser;
  static double cof[6]={76.18009172947146,-86.50532032941677,24.01409824083091,-1.231739572450155,0.1208650973866179e-2,-0.5395239384953e-5};
  int j;
  y = x = xx;
  tmp = x+5.5;
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