#pragma once

class EmpiricalSelACProfileBidimArray : public SimpleBidimArray<vector<double> > {

    public:
    EmpiricalSelACProfileBidimArray(const vector<vector<double>>& inaafitness, const vector<double>& inG, double inpsi) :
        SimpleBidimArray<vector<double> > (Naa, inG.size(), vector<double>(Naa,1.0/Naa)),
        aafitness(inaafitness), G(inG), psi(inpsi) {
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
            double tmp = psi * G[j] * log(aafitness[i][a]);
            v[a] = 1e-50;
            if (tmp < 100)  {
                v[a] += exp(tmp);
            }
            tot += v[a];
        }
        for (int a=0; a<Naa; a++)   {
            v[a] /= tot;
        }
    }

  private:

    const vector<vector<double>>& aafitness;
    const vector<double>& G;
    double psi;
};

class UnconstrainedSelACProfileBidimArray : public SimpleBidimArray<vector<double> > {

    public:
    UnconstrainedSelACProfileBidimArray(const Selector<vector<double>>& inaafitness, const vector<double>& inG, double inpsi = 1) :
        SimpleBidimArray<vector<double> > (inaafitness.GetSize(), inG.size(), vector<double>(Naa,1.0/Naa)),
        aafitness(inaafitness), G(inG), psi(inpsi) {
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
            double tmp = psi * G[j] * log(aafitness.GetVal(i)[a]);
            v[a] = 1e-50;
            if (tmp < 100)  {
                v[a] += exp(tmp);
            }
            tot += v[a];
        }
        for (int a=0; a<Naa; a++)   {
            v[a] /= tot;
        }
    }

  private:

    const Selector<vector<double>>& aafitness;
    const vector<double>& G;
    double psi;
};

