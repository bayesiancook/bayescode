#pragma once

class BidimSum {
    public:
        BidimSum(const BidimSelector<double>& inx1, const BidimSelector<double>& inx2):
            x1(inx1), x2(inx2) {}

        ~BidimSum() {}

        int GetNrow() const {return x1.GetNrow();}
        int GetNcol() const {return x1.GetNcol();}

        double GetVal(int i, int j) const {
            return x1.GetVal(i,j) + x2.GetVal(i,j);
        }

    private:
        const BidimSelector<double>& x1;
        const BidimSelector<double>& x2;
};

