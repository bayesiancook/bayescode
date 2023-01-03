#pragma once

class BidimWeightedSum {
    public:
        BidimWeightedSum(const BidimSelector<double>& inx1, 
                const BidimSelector<double>& inx2, 
                const Selector<double>& intimescale, 
                double& inmean2,
                int inrelative) : 
            x1(inx1), x2(inx2), timescale(intimescale), mean2(inmean2), relative(inrelative) {}

        ~BidimWeightedSum() {}

        int GetNrow() const {return x1.GetNrow();}
        int GetNcol() const {return x1.GetNcol();}

        double GetVal(int i, int j) const {
            if (relative)   {
                return x1.GetVal(i,j) + mean2 * x2.GetVal(i,j) / timescale.GetVal(j);
            }
            return x1.GetVal(i,j) + mean2 * x2.GetVal(i,j);
        }

    private:
        const BidimSelector<double>& x1;
        const BidimSelector<double>& x2;
        const Selector<double>& timescale;
        double& mean2;
        int relative;
};

