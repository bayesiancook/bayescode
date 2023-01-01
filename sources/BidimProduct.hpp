
#pragma once

class BidimProduct {
    public:
        BidimProduct(const Selector<double>& in_rows, const Selector<double>& in_cols) : 
            rows(in_rows), cols(in_cols) {}

        ~BidimProduct() {}

        int GetNrow() const {return rows.GetSize();}
        int GetNcol() const {return cols.GetSize();}

        double GetVal(int i, int j) const {return rows.GetVal(i) * cols.GetVal(j);}

    private:
        const Selector<double>& rows;
        const Selector<double>& cols;
};

