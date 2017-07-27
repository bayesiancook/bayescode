#ifndef DISCBETA_H
#define DISCBETA_H

#include "Array.hpp"
#include "cdf.hpp"

class DiscBetaWithPos : public SimpleArray<double>  {

    public:

    DiscBetaWithPos(int inncat, double inalpha, double inbeta, double inposw, double inposom) : SimpleArray<double>(inncat+1) , ncat(inncat), weight(inncat+1), alpha(inalpha), beta(inbeta), posw(inposw), posom(inposom) {
            ComputeDiscBeta();
            (*this)[ncat] = posom;
            ComputeWeights();
    }

    ~DiscBetaWithPos() {}

    void SetAlphaBeta(double inalpha, double inbeta)    {
        alpha = inalpha;
        beta = inbeta;
        ComputeDiscBeta();
    }
    
    void SetPos(double inposw, double inposom)  {
        posw = inposw;
        posom = inposom;
        (*this)[ncat] = posom;
        ComputeWeights();
    }

    const vector<double>& GetWeights() const {return weight;}

    double GetPostProbArray(const OmegaSuffStat& suffstat, vector<double>& postprob) const {

            double logp[ncat+1];
            double max = 0;
            for (int cat=0; cat<=ncat; cat++)  {
                logp[cat] = suffstat.GetLogProb(GetVal(cat));
                if ((!cat) || (max < logp[cat]))    {
                    max = logp[cat];
                }
            }
            double tot = 0;
            for (int cat=0; cat<=ncat; cat++)  {
                postprob[cat] = weight[cat] * exp(logp[cat] - max);
                tot += postprob[cat];
            }
            for (int cat=0; cat<=ncat; cat++)  {
                postprob[cat] /= tot;
            }
            double ret = log(tot) + max;
            if (isinf(ret)) {
                cerr << "ret is inf: " << tot << '\t' << max << '\n';
                exit(1);
            }
            if (isnan(ret)) {
                cerr << "ret is nan: " << tot << '\t' << max << '\n';
                for (int cat=0; cat<=ncat; cat++)   {
                    cerr << GetVal(cat) << '\t' << weight[cat] << '\n';
                }
                exit(1);
            }
            return ret;
    }

    double GetPostProbArray(OmegaSuffStatArray& suffstatarray, vector<vector<double> >& postprobarray)  const {

        double total = 0;
        for (int i=0; i<suffstatarray.GetSize(); i++)   {
            total += GetPostProbArray(suffstatarray.GetVal(i),postprobarray[i]);
        }
        return total;
    }

    private:

    // still numerically unstable
    void ComputeDiscBeta()   {
        for (int cat=0; cat<ncat; cat++)  {
            (*this)[cat] = invbetaInc(alpha,beta,((double) (cat+0.5))/ncat);
        }
    }

    void ComputeWeights()   {
        if (posw > 1)   {
            cerr << "error: pos weight greater than 1\n";
            exit(1);
        }
        if (posw < 0)   {
            cerr << "error: negative pos weight\n";
            exit(1);
        }
        for (int cat=0; cat<ncat; cat++)    {
            weight[cat] = (1-posw)/ncat;
        }
        weight[ncat] = posw;
    }

    int ncat;
    vector<double> weight;
    double alpha;
    double beta;
    double posw;
    double posom;
};

#endif

