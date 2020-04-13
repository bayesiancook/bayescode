
#pragma once

class SelACProfileBidimArray : public SimpleBidimArray<vector<double> > {

    public:
    SelACProfileBidimArray(const vector<double>& inaadist, const vector<double>& inG, double inpsi) :
        SimpleBidimArray<vector<double> > (Naa, inG.size(), vector<double>(Naa,1.0/Naa)),
        aadist(inaadist), G(inG), psi(inpsi) {
            Update();
    }

    void SetPsi(double inpsi)   {
        psi = inpsi;
    }

    void Update()   {
        for (int i = 0; i < GetNrow(); i++) {
            for (int j = 0; j < GetNcol(); j++) {
                Update(i,j);
            }
        }
    }

    void UpdateRow(int i)   {
        if (i > GetNrow())  {
            cerr << "error in selac array: row overflow\n";
            exit(1);
        }
        for (int j = 0; j < GetNcol(); j++) {
            Update(i,j);
        }
    }

    void UpdateColumn(int j)   {
        if (j > GetNcol())  {
            cerr << "error in selac array: col overflow\n";
            exit(1);
        }
        for (int i = 0; i < GetNrow(); i++) {
            Update(i,j);
        }
    }

    void Update(int i, int j)   {
        vector<double>& v = (*this)(i,j);
        double tot = 0;
        for (int a=0; a<Naa; a++)   {
            if (a == i) {
                v[a] = 1.0;
            }
            else    {
                double tmp = psi * G[j] * aadist[rrindex(a,i)];
                v[a] = 1e-8;
                if (tmp < 100)  {
                    v[a] += exp(-psi*G[j]*aadist[rrindex(a,i)]);
                }
            }
            tot += v[a];
        }
        for (int a=0; a<Naa; a++)   {
            v[a] /= tot;
        }
    }

  private:

    const vector<double>& aadist;
    const vector<double>& G;
    double psi;

    int rrindex(int i, int j) const {
        return (i < j) ? (2 * Naa - i - 1) * i / 2 + j - i - 1
                       : (2 * Naa - j - 1) * j / 2 + i - j - 1;
    }
};

