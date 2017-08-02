#ifndef DISCBETA_H
#define DISCBETA_H

#include "Array.hpp"
#include "cdf.hpp"

class DiscBetaWithPos : public SimpleArray<double>  {

    public:

    DiscBetaWithPos(int inncat, double inalpha, double inbeta, double inposw, double inposom, vector<double>& inw) : SimpleArray<double>(inncat+3) , ncat(inncat), weight(inncat+3), alpha(inalpha), beta(inbeta), posw(inposw), posom(inposom), w(inw) {
            ComputeDiscBeta();
            (*this)[0] = 1e-2;
            (*this)[ncat+1] = 1;
            (*this)[ncat+2] = posom;
            ComputeWeights();
    }

    ~DiscBetaWithPos() {}

    void SetAlphaBeta(double inalpha, double inbeta)    {
        alpha = inalpha;
        beta = inbeta;
        ComputeDiscBeta();
    }
    
    void SetPurifWeight(const vector<double>& inw)    {
        w = inw;
        ComputeWeights();
    }

    void SetPos(double inposw, double inposom)  {
        posw = inposw;
        posom = inposom;
        (*this)[ncat+2] = posom;
        ComputeWeights();
    }

    const vector<double>& GetWeights() const {return weight;}

    double GetPostProbArray(const OmegaSuffStat& suffstat, vector<double>& postprob) const {

            double logp[GetSize()];
            double max = 0;
            for (int cat=0; cat<GetSize(); cat++)  {
                logp[cat] = suffstat.GetLogProb(GetVal(cat));
                if ((!cat) || (max < logp[cat]))    {
                    max = logp[cat];
                }
            }
            double tot = 0;
            for (int cat=0; cat<GetSize(); cat++)  {
                postprob[cat] = weight[cat] * exp(logp[cat] - max);
                tot += postprob[cat];
            }
            for (int cat=0; cat<GetSize(); cat++)  {
                postprob[cat] /= tot;
            }
            double ret = log(tot) + max;
            if (isinf(ret)) {
                cerr << "ret is inf: " << tot << '\t' << max << '\n';
                exit(1);
            }
            if (isnan(ret)) {
                cerr << "ret is nan: " << tot << '\t' << max << '\n';
                for (int cat=0; cat<GetSize(); cat++)   {
                    cerr << GetVal(cat) << '\t' << weight[cat] << '\t' << logp[cat] << '\t' << postprob[cat] << '\n';
                }
                cerr << tot << '\t' << log(tot) << '\t' << max << '\n';
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
        if (ncat == 1)  {
            (*this)[1] = alpha / (alpha+beta);
        }
        else    {
            for (int cat=0; cat<ncat; cat++)  {
                (*this)[cat+1] = invbetaInc(alpha,beta,((double) (cat+0.5))/ncat);
            }
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
            weight[cat+1] = w[1]*(1 - posw) / ncat;
        }
        weight[0] = w[0]*(1 - posw);
        weight[ncat+1] = w[2]*(1 - posw);
        weight[ncat+2] = posw;
    }

    int ncat;
    vector<double> weight;
    double alpha;
    double beta;
    double posw;
    double posom;
    vector<double> w;
};

#endif

