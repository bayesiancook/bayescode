#ifndef INCBETA_H
#define INCBETA_H

// numerical recipe for beta function
// taken from statweb.standford.edu/~serban/gxna/src/cdf.cpp
// (Nacu et al, Bioinformatics 2007, 25:850)
// inverse beta function: dichotomy

#include "Random.hpp"

double logGamma(double x);
double betaContFrac(double a, double b, double x);
double betaInc(double a, double b, double x);
double invbetaInc(double a, double b, double p);
/*
double normCDF(const double x);
double normCDFInv(const double x);
double tCDF(const double x, const double n);
double fCDF(const double x, const double n1, const double n2);
double zfCDF(const double x, const double n1, const double n2);
double ztCDF(const double x, const double n);
double digamma(const double x);
double trigamma(const double x);
double trigammainv(const double y);
*/

#endif
