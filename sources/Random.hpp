#ifndef RANDOM_H
#define RANDOM_H

#include "Eigen/Dense"
//
// c++11
// using EMatrix = Eigen::MatrixXd;
// using EVector = Eigen::VectorXd;

typedef Eigen::MatrixXd EMatrix;
typedef Eigen::VectorXd EVector;

using namespace std;

#define MT_LEN 624  // (VL) required for magic
#include <vector>

const double Pi = 3.1415926535897932384626;

class Random {
  public:
    static const double INFPROB;

    Random(int seed = -1);

    static void InitRandom(int seed = -1);

    static int GetSeed();

    static double Uniform();
    static int ApproxBinomial(int N, double p);
    static int Poisson(double mu);
    static double Gamma(double alpha, double beta);
    static double sNormal();
    static double sExpo();
    static double sGamma(double);
    static double sGammanew(double);

    static int Choose(int);
    static int FiniteDiscrete(int n, const double* probarray);
    static void DrawFromUrn(int*, int n, int N);
    static int DrawFromDiscreteDistribution(const EVector& prob, int nstate);
    static int DrawFromDiscreteDistribution(const double* prob, int nstate);
    static int DrawFromDiscreteDistribution(const std::vector<double>& prob);

    static double logGamma(double alpha);

    static double logMultivariateGamma(double a, int p);

    static double ProfileProposeMove(std::vector<double>& profile, int dim, double tuning, int n);
    static double RealVectorProposeMove(std::vector<double>& x, int dim, double tuning, int n);

    static double GetEntropy(const std::vector<double>& profile);

    static double GammaSample(double alpha, double beta);
    static double BetaSample(double alpha, double beta);
    static void DirichletSample(vector<double>& x, const vector<double>& center, double concentration);

    static double logGammaDensity(double x, double alpha, double beta);
    static double logBetaDensity(double x, double alpha, double beta);
    static double logDirichletDensity(const vector<double>& x, const vector<double>& center, double concentration);

  private:
    static int Seed;
    static int mt_index;
    static unsigned long long mt_buffer[MT_LEN];
};

#endif  // RANDOM_H
