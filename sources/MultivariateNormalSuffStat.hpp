#pragma once

#include "SuffStat.hpp"
#include "CovMatrix.hpp"

class MultivariateNormalSuffStat : public SuffStat  {

    public:

    MultivariateNormalSuffStat(int indim) : covmat(indim), n(0) {}

    int GetDim() const  {
        return covmat.GetDim();
    }

    void Reset()    {
        for (int i=0; i<GetDim(); i++)  {
            for (int j=0; j<GetDim(); j++)  {
                covmat.setval(i,j,0);
            }
        }
        n = 0;
    }

    CovMatrix covmat;
    int n;
};
