
#include "Random.hpp"
#include <iostream>

double logscore(double lambda, double mu)   {

    int n = 2;
    double tot = (n-1)*log(2.0 * lambda) + n*log((lambda-mu)*(lambda-mu));
    /*
    double temp = lambda-mu*exp(-(lambda-mu));
    double f = -(lambda-mu) - log(temp*temp);
    tot += f;
    */

    double temp0 = lambda-mu*exp(-2*(lambda-mu));
    double f0 = -2*(lambda-mu) - log(temp0*temp0);
    tot += f0;

    double temp1 = lambda-mu*exp(-(lambda-mu));
    double f1 = -(lambda-mu) - log(temp1*temp1);
    tot += f1;

    tot -= lambda + mu;
    return tot;
}


int main()  {

    double lambda = 1.0;
    double mu = 1.0;

    int nrep = 1000000;
    double delta = 1;
    int burnin = 1000;
    int nsample = 0;
    double meanlambda = 0;
    double meanmu = 0;

    for (int rep=0; rep<nrep; rep++)    {

            double bklambda = lambda;
            double bkmu = mu;

            double logprob1 = logscore(lambda,mu);

            double m = delta * (Random::Uniform() - 0.5);
            double e = exp(m);

            if (Random::Uniform() < 0.5)    {
                lambda *= e;
            }
            else    {
                mu *= e;
            }

            double logprob2 = logscore(lambda,mu);
            double logratio = logprob2 - logprob1 + m;

            if (log(Random::Uniform()) > logratio)  {
                lambda = bklambda;
                mu = bkmu;
            }


            if (rep > burnin)   {
                meanlambda += lambda;
                meanmu += mu;
                nsample++;
            }
    }

    meanlambda /= nsample;
    meanmu /= nsample;

    std::cout << meanlambda << '\t' << meanmu << '\n';

}


