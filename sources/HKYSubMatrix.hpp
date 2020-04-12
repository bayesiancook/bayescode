#pragma once

#include "SubMatrix.hpp"

class HKYSubMatrix : public virtual SubMatrix   {

  public:

    HKYSubMatrix(double inkappa, double ingc, bool innormalise = false) : 
        SubMatrix(4,innormalise), kappa(inkappa), gc(ingc)  {
    }

    ~HKYSubMatrix() override = default;

    void SetKappa(double inkappa)   {
        kappa = inkappa;
        CorruptMatrix();
    }

    void SetGC(double ingc) {
        gc = ingc;
        CorruptMatrix();
    }

  protected:

    void ComputeStationary() const override {
        mStationary[0] = (1-gc) / 2;
        mStationary[1] = gc/2;
        mStationary[2] = gc/2;
        mStationary[3] = (1-gc) / 2;
    }

    void ComputeArray(int i) const override {
        double total = 0;
        for (int j=0; j<Nstate; j++)	{
            if (i!=j)	{
                Q(i,j) = Stationary(j);
                if (((i==0) && (j==2)) || ((i==2) && (j==0)) || ((i==1) && (j==3)) || ((i==3) && (j==1)))	{
                    Q(i,j) *= kappa;
                }
                total += Q(i,j);
            }
        }
        Q(i,i) = -total;
    }

  private:
    double kappa;
    double gc;
};

