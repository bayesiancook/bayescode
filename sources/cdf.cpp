
#include "cdf.hpp"
#include <math.h>
#include <iostream>
#include <cstdlib>
using namespace std;

// const double PI = 3.1415926535;
const double SQRT2PI = 2.5066282746;
const double SQRT2 = 1.414214;

// log gamma using the Lanczos approximation
// must have x > 0
double logGamma(double x) {
  const double c[8] = { 676.5203681218851, -1259.1392167224028,
			 771.32342877765313, -176.61502916214059, 
			 12.507343278686905, -0.13857109526572012, 
			 9.9843695780195716e-6, 1.5056327351493116e-7 };
  double sum = 0.99999999999980993;
  double y = x;
  for (int j = 0; j < 8; j++)
    sum += c[j] / ++y;
  return log(SQRT2PI * sum / x) - (x + 7.5) + (x + 0.5) * log(x + 7.5);
}

// helper function for incomplete beta
// computes continued fraction
// source: Numerical Recipes in C
double betaContFrac(double a, double b, double x) {
  const int MAXIT = 1000;
  const double EPS = 3e-7;
  const double FPMIN = 1e-30;
  double qab = a + b;
  double qap = a + 1;
  double qam = a - 1;
  double c = 1;
  double d = 1 - qab * x / qap;
  if (fabs(d) < FPMIN) d = FPMIN;
  d = 1 / d;
  double h = d;
  int m;
  for (m = 1; m <= MAXIT; m++) {
    int m2 = 2 * m;
    double aa = m * (b-m) * x / ((qam + m2) * (a + m2));
    d = 1 + aa * d;
    if (fabs(d) < FPMIN) d = FPMIN;
    c = 1 + aa / c;
    if (fabs(c) < FPMIN) c = FPMIN;
    d = 1 / d;
    h *= (d * c);
    aa = -(a+m) * (qab+m) * x / ((a+m2) * (qap+m2));
    d = 1 + aa * d;
    if (fabs(d) < FPMIN) d = FPMIN;
    c = 1 + aa / c;
    if (fabs(c) < FPMIN) c = FPMIN;
    d = 1 / d;
    double del = d*c;
    h *= del;
    if (fabs(del - 1) < EPS) break;
  }
  if (m > MAXIT) {
    cerr << "betaContFrac: too many iterations\n";
  }
  return h;
}

// incomplete beta function
// must have 0 <= x <= 1
double betaInc(double a, double b, double x) {
  if (x == 0)
    return 0;
  else if (x == 1)
    return 1;
  else {
    double logBeta = logGamma(a+b) - logGamma(a) - logGamma(b)
      + a * log(x) + b * log(1-x);
    if (x < (a+1) / (a+b+2))
      return exp(logBeta) * betaContFrac(a, b, x) / a;
    else
      return 1 - exp(logBeta) * betaContFrac(b, a, 1-x) / b;
  }
}

double invbetaInc(double a, double b, double p) {

    double x1 = 0;
    double x2 = 1;
    double eps = 1e-10;
    while (x2-x1 > eps) {
        double xmid = (x1+x2)/2;
        double pmid = betaInc(a,b,xmid);
        if (p < pmid)   {
            x2 = xmid;
        }
        else    {
            x1 = xmid;
        }
    }
    return (x1+x2)/2;
}

/*
double invbetaInc(double a, double b, double p) {

    double a1 = a-1;
    double b1 = b-1;
    double error = 1e-8;

    if (p<=0)    {
       return 0;
    }

    if (p>=1)    {
       return 1;
    }

    double x = 0;

    if ((a >= 1.0) && (b >= 1.0))   {
        double pp = 0;
        if (p<0.5)  {
            pp = p;
        }
        else    {
            pp = 1-p;
        }
        double t = sqrt(-2*log(pp));
        x = (2.30753 + t*0.27061)/(1.0 + t * (0.99229 + t*0.04481)) - t;
        if (p < 0.5)    {
            x = -x;
        }
        a1 = (x*x-3)/6;
        double h = 2.0/(1.0/(2.0*a-1.0) + 1.0/(2.0*b-1.0));
        double w = (x*sqrt(a1+h)/h) - (1.0/(2.0*b-1.0) - 1.0/(2.0*a-1.0)) * (a1 + 5.0/6.0 - 2.0/(3.0*h));
        x = a/(a + b*exp(2*w));
    }
    else    {

        double lna = log(a / (a+b));
        double lnb = log(b / (a+b));
        double t = exp(a*lna)/a;
        double u = exp(b*lnb)/b;
        double w = t+u;
        if (p < t/w)    {
            x = exp(1.0/a*log(a*w*p));
        }
        else    {
            x = 1.0 - exp(1.0/b*log(b*w*(1-p)));
        }
    }

    double afac = -logGamma(a) - logGamma(b) + logGamma(a+b);
    int nrep = 1000;
    int i = 0;
    double t = 0;
    do  {
        i++;
        if ((x == 0) || (x == 1.0)) {
            return x;
        }
        double err = betaInc(a,b,x) - p;
        t = exp(a1*log(x) + b1*log(1-x) + afac);
        double u = err/t;
        double min = 1;
        double tmp = u*(a1/x - b1/(1-x));
        if (min > tmp)  {
            min = tmp;
        }

        t = u / (1.0 - 0.5*min);
        x-=t;
        if (x <= 0) {
            x = 0.5*(x+t);
        }
        if (x >= 1.0)   {
            x = 0.5*(x+t+1);
        }
    }
    while ((i < nrep) && (fabs(t) < error*x));

    return x;
}

// error integral using the Bagby approximation
// max error has order 10^-5, should replace it by something better
double normCDF(const double x) {
  double z1 = exp(-x*x/2);
  double z2 = exp(-x*x*0.5857864); // const is 2-sqrt(2)
  double s = 7*z1 + 16*z2 + (7 + 0.7853982*x*x)*z1*z1; // const is pi/4
  double t = sqrt(1 - s/30) / 2;
  return x > 0 ? t + 0.5 : 0.5 - t;
}

// inverse error function
// gives quantiles for gaussian: f(0)=-Inf, f(1/2)=0, f(1)=Inf
// overflows within 10^-10 or so of 0 and 1
// max error has order 10^-3, should replace it by something better
double normCDFInv(const double x) {
  const double a = 0.1400123;
  const double b = 4.546885; // 2 / (pi * a)
  double ln = 0.5 * (log(4.0) + log(x) + log(1-x));
  double erfi = sqrt(-b - ln + sqrt((b+ln)*(b+ln) - 2*ln/a));
  return x >= 0.5 ? erfi * SQRT2 : -erfi * SQRT2;
}

double tCDF(const double x, const double n) {
  return 1 - 0.5 * betaInc(n / 2, 0.5, n / (n + x * x));
}

double fCDF(const double x, const double n1, const double n2) {
  return 1 - betaInc(n2 / 2, n1 / 2, n2 / (n2 + n1 * x));
}

// normal equivalent to f statistic, avoids overflow
// returns y with P(|normal| < y) = P(F < x)
// needs x >= 0
double zfCDF(const double x, const double n1, const double n2) {
  double y = n2 / (n2 + n1 * x);
  double a = n2 / 2;
  double b = n1 / 2;
  double bInc, logBInc;
  if (y >= (a+1) / (a+b+2)) {
    bInc = betaInc(a, b, y);
    logBInc = log(bInc);
  }
  else {
    double logBeta = logGamma(a+b) - logGamma(a) - logGamma(b)
      + a * log(y) + b * log(1-y);
    logBInc = logBeta + log(betaContFrac(a, b, y)) - log(a);
    bInc = exp(logBInc);
  }
  const double A = 0.1400123;
  const double B = 4.546885; // 2 / (pi * A)
  double ln = 0.5 * (logBInc + log(2-bInc));
  double erfi = sqrt(-B - ln + sqrt((B+ln)*(B+ln) - 2*ln/A));
  return erfi * SQRT2;
}

// normal equivalent to t statistic, avoids overflow
// returns y with P(|normal| < y) = P(T < x)
// needs x >= 0
double ztCDF(const double x, const double n) {
  return x >= 0 ? zfCDF(x * x, 1, n) : -zfCDF(x * x, 1, n);
}

// digamma and trigamma have slow and sloppy but easy implementations
double digamma(const double x) {
  if (x < 20)
    return digamma(x+1) - 1/x;
  else if (x > 21)
    return digamma(x-1) + 1/(x-1);
  else // interpolate
    return 2.970524 * (21-x) + 3.020524 * (x-20);
}

double trigamma(const double x) {
  if (x < 20)
    return trigamma(x+1) + 1/x/x;
  else if (x > 21)
    return trigamma(x-1) - 1/(x-1)/(x-1);
  else // interpolate
    return 0.05127082 * (21-x) + 0.04877082 * (x-20);
}

double trigammainv(const double y) {
  double x1 = 1 / y + 0.5;
  double x2 = 1 / (y + 1 / x1) - 0.5;
  return (x1 + x2) / 2;
}
*/
