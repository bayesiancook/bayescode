#include "Random.hpp"
#include "grantham.hpp"
#include "LG.hpp"
#include "MeanSubMatrix.hpp"
#include "HKYSubMatrix.hpp"
#include "SelACProfileBidimArray.hpp"
#include "AAMutSelPolyMutCodonMatrixBidimArray.hpp"
#include "MeanSubMatrix.hpp"
#include "AAProjectedCodonSubMatrix.hpp"
#include "ProbModel.hpp"
#include <fstream>

class FastEmpiricalSelACModel : public ProbModel {

    public:

    string name;

    double epsilon;

    double kappa;
    double gamma;
    double u;
    HKYSubMatrix *nucmatrix;

    CodonStateSpace* statespace;

    // amino acid selection
    int Gcat;
    double Ginvshape;
    vector<double> G;
    vector<double> Gweight;
    double psi;

    vector<double> obsweight;
    vector<vector<double>> obscompaafreq;
    vector<vector<double>> predcompaafreq;

    vector<double> aadist;
    vector<double> granthaadist;
    int Naa;
    int Naarr;

    // putting things together
    SelACProfileBidimArray* selacprofiles;

    // to do with double and triple mutants
    AAMutSelPolyMutCodonMatrixBidimArray* codonmatrices;
    MeanCodonSubMatrixFromBidimArray* meancodonmatrix;
    AAProjectedCodonSubMatrix* aamatrix;

    // mean dN/dS across sites
    double obsmeanomega;
    double predmeanomega;

    // variance in dN/dS across sites
    double obsrelvaromega;
    double predrelvaromega;

    vector<double> aadist_acc;
    vector<double> aadist_tot;
    vector<double> G_acc;
    vector<double> G_tot;
    vector<double> psi_acc;
    vector<double> psi_tot;
    vector<double> gpsi_acc;
    vector<double> gpsi_tot;

    vector<vector<double>> preddnds;

    // allocation

    FastEmpiricalSelACModel(string aafitnessfile, int inGcat, double inkappa, double ingamma, double inu, double inobsmeanomega, double inobsrelvaromega, double inepsilon, string inname) {

        name = inname;
        Gcat = inGcat;
        obsmeanomega = inobsmeanomega;
        obsrelvaromega = inobsrelvaromega;
        epsilon = inepsilon;

        kappa = inkappa;
        gamma = ingamma;
        u = inu;

        Naa = 20;
        Naarr = Naa * (Naa-1) / 2;

        preddnds.assign(Naa, vector<double>(Gcat,0));

        obsweight.assign(Naa,0);
        obscompaafreq.assign(Naa, vector<double>(Naa,0));
        predcompaafreq.assign(Naa, vector<double>(Naa,0));
        ifstream is(aafitnessfile.c_str());
        int naa;
        is >> naa;
        if (naa != Naa) {
            cerr << "error when reading aafitness file\n";
            exit(1);
        }
        for (int i=0; i<Naa; i++)   {
            is >> obsweight[i];
            for (int j=0; j<Naa; j++)   {
                is >> obscompaafreq[i][j];
                if (! obscompaafreq[i][j])  {
                    obscompaafreq[i][j] = 1e-7;
                }
            }
        }

        // monitoring
        aadist_acc.assign(3,0);
        aadist_tot.assign(3,0.1);
        G_acc.assign(3,0);
        G_tot.assign(3,0.1);
        psi_acc.assign(3,0);
        psi_tot.assign(3,0.1);
        gpsi_acc.assign(4,0);
        gpsi_tot.assign(4,0.1);

        nucmatrix = new HKYSubMatrix(kappa,gamma);

        Ginvshape = 2.0;
        G.assign(Gcat, 1.0/Gcat);
        UpdateG();
        psi = 3.5;
        statespace = new CodonStateSpace(Universal);

        aadist.assign(Naarr, 1.0);
        UpdateGrantham(grantham_wcom_selac, grantham_wpol, grantham_wvol);
        granthaadist.assign(Naarr, 1.0);
        for (int i=0; i<Naarr; i++) {
            granthaadist[i] = aadist[i];
        }
        /*
        int index = 0;
        for (int i=0; i<Naa; i++)   {
            for (int j=i+1; j<Naa; j++)   {
                aadist[index] = -log(obscompaafreq[i][j] + obscompaafreq[j][i] + 0.01);
                cerr << AminoAcids[i] << AminoAcids[j] << '\t' << aadist[index] << '\n';
                index++;
            }
        }
        */

        Gweight.assign(Gcat, 1.0/Gcat);

        selacprofiles = new SelACProfileBidimArray(aadist,G,psi);

        codonmatrices = new AAMutSelPolyMutCodonMatrixBidimArray(*selacprofiles, *statespace, *nucmatrix, u, 1.0);

        Update();
    }

    // updates and logprob 

    void UpdateG()  {
        Random::DiscGamma(G,1.0/Ginvshape);
    }

    int rrindex(int i, int j) const {
        return (i < j) ? (2 * Naa - i - 1) * i / 2 + j - i - 1
                       : (2 * Naa - j - 1) * j / 2 + i - j - 1;
    }

    void UpdateGrantham(double wcom, double wpol, double wvol)   {
        double tot = 0;
        for (int a=0; a<Naa; a++)   {
            for (int b=a+1; b<Naa; b++)   {
                double tcom = grantham_com[b] - grantham_com[a];
                double tpol = grantham_pol[b] - grantham_pol[a];
                double tvol = grantham_vol[b] - grantham_vol[a];
                double d = sqrt(wcom*tcom*tcom + wpol*tpol*tpol + wvol*tvol*tvol);
                aadist[rrindex(a,b)] = d;
                tot += d;
            }
        }
        tot /= Naarr;
        for (int a=0; a<Naa; a++)   {
            for (int b=a+1; b<Naa; b++)   {
                aadist[rrindex(a,b)] /= tot;
            }
        }
    }

    void GetPredictedDNDS(double& mean, double& relvar) const  {
        double m0 = 0;
        double m1 = 0;
        double m2 = 0;
        for (int i=0; i<Naa; i++) {
            for (int j=0; j<Gcat; j++)  {
                double dnds = preddnds[i][j];
                m0 += obsweight[i] * Gweight[j];
                m1 += obsweight[i] * Gweight[j] * dnds;
                m2 += obsweight[i] * Gweight[j] * dnds * dnds;
            }
        }
        m1 /= m0;
        m2 /= m0;
        m2 -= m1*m1;
        mean = m1;
        relvar = m2 / m1 / m1;
    }


    void ComputePredictedDNDS() {
        for (int a=0; a<Naa; a++)   {
            ComputePredictedDNDS(a);
        }
    }

    void ComputePredictedDNDS(int a)    {
        for (int j=0; j<Gcat; j++)  {
            preddnds[a][j] = codonmatrices->GetVal(a,j).GetPredictedDNDS();
        }
    }

    void ComputeSummaryStatistics() {
        GetPredictedDNDS(predmeanomega, predrelvaromega);
        codonmatrices->GetMeanAAFrequencies(predcompaafreq);
    }

    void Update() override {
        UpdateG();
        selacprofiles->SetPsi(psi);
        selacprofiles->Update();
        codonmatrices->Corrupt();
        ComputePredictedDNDS();
        ComputeSummaryStatistics();
    }

    void Update(int a, int b)  {
        selacprofiles->UpdateRow(a);
        selacprofiles->UpdateRow(b);
        codonmatrices->CorruptRow(a);
        codonmatrices->CorruptRow(b);
        ComputePredictedDNDS(a);
        ComputePredictedDNDS(b);
        ComputeSummaryStatistics();
    }

    double GetMSD() const {
        double diff = 0;
        double tmp = 50*(obsmeanomega - predmeanomega)/obsmeanomega;
        diff += tmp*tmp;
        tmp = 50*(obsrelvaromega - predrelvaromega)/obsrelvaromega;
        diff += tmp*tmp;

        for (int i=0; i<Naa; i++)   {
            for (int j=0; j<Naa; j++)   {
                double tmp = log(obscompaafreq[i][j]) - log(predcompaafreq[i][j]);
                diff += tmp*tmp;
            }
        }
        return diff;
    }

    double AALogPrior() const   {
        double tot = 0;
        for (int i=0; i<Naarr; i++) {
            tot -= aadist[i];
        }
        return tot;
    }

    double GetLogPrior() const {
        double tot = 0;
        tot -= psi;
        tot -= Ginvshape;
        tot += AALogPrior();
        return tot;
    }

    double GetLogProb() const override {
        double diff = GetMSD();
        return GetLogPrior() - diff * diff / 2 / epsilon / epsilon;
    }

    // moving

    double Move() override {
        int nrep = 10;
        for (int rep=0; rep<nrep; rep++)    {

            aadist_acc[0] += FastMoveAADist(100,1.0);
            aadist_tot[0] ++;
            aadist_acc[1] += FastMoveAADist(100,0.1);
            aadist_tot[1] ++;

            MoveGinvshape();
            MovePsi();

            /*
            gpsi_acc[0] += GPsiCompMove(3, 1.0, 40.0);
            gpsi_tot[0] ++;
            gpsi_acc[1] += GPsiCompMove(3, 0.1, 40.0);
            gpsi_tot[1] ++;
            gpsi_acc[2] += GPsiCompMove(3, 1.0, 10.0);
            gpsi_tot[2] ++;
            gpsi_acc[3] += GPsiCompMove(3, 0.1, 10.0);
            gpsi_tot[3] ++;
            */

            AAPsiCompMove(1,10);
            AAPsiCompMove(0.1,10);
        }
        return 1.0;
    }

    double FastMoveAADist(int nrep, double tuning) {
        double nacc = 0;
        double ntot = 0;
        for (int i=0; i<nrep; i++)  {
            int a = (int) (Naa * Random::Uniform());
            int b = (int) ((Naa-1) * Random::Uniform());
            if (b >= a) b++;
            int j = rrindex(a,b);
            double delta = -GetLogProb();
            double m = tuning * (Random::Uniform() - 0.5);
            double e = exp(m);
            aadist[j] *= e;
            Update(a,b);
            delta += GetLogProb();
            delta += m;
            int acc = (log(Random::Uniform()) < delta);
            if (acc)    {
                nacc++;
            }
            else    {
                aadist[j] /= e;
                Update(a,b);
            }
            ntot++;
        }
        return nacc / ntot;
    }

    double GPsiCompMove(int nrep, double tuning, double f)   {

        double nacc = 0;
        double ntot = 0;
        for (int rep=0; rep<nrep; rep++)    {
            double delta = -GetLogProb();
            double m = tuning * (Random::Uniform() - 0.5);
            psi += m;
            Ginvshape += m/f;
            if ((psi > 0) && (Ginvshape > 0))   {
                Update();
                delta += GetLogProb();
                int acc = (log(Random::Uniform()) < delta);
                if (acc)    {
                    nacc++;
                }
                else    {
                    psi -= m;
                    Ginvshape -= m/f;
                    Update();
                }
            }
            else    {
                psi -= m;
                Ginvshape -= m/f;
            }
            ntot++;
        }
        return nacc/ntot;
    }

    double MoveGinvshape()  {
        G_acc[0] += ScalingMove(Ginvshape, 0.3, 3, &FastEmpiricalSelACModel::GetLogProb, &FastEmpiricalSelACModel::Update, this);
        G_tot[0] ++;
        G_acc[1] += ScalingMove(Ginvshape, 0.03, 3, &FastEmpiricalSelACModel::GetLogProb, &FastEmpiricalSelACModel::Update, this);
        G_tot[1] ++;
        G_acc[2] += ScalingMove(Ginvshape, 0.003, 3, &FastEmpiricalSelACModel::GetLogProb, &FastEmpiricalSelACModel::Update, this);
        G_tot[2] ++;
        return 1.0;
    }

    double MovePsi()  {
        psi_acc[0] += ScalingMove(psi, 0.3, 3, &FastEmpiricalSelACModel::GetLogProb, &FastEmpiricalSelACModel::Update, this);
        psi_tot[0] ++;
        psi_acc[1] += ScalingMove(psi, 0.03, 3, &FastEmpiricalSelACModel::GetLogProb, &FastEmpiricalSelACModel::Update, this);
        psi_tot[1] ++;
        psi_acc[2] += ScalingMove(psi, 0.003, 3, &FastEmpiricalSelACModel::GetLogProb, &FastEmpiricalSelACModel::Update, this);
        psi_tot[2] ++;
        return 1.0;
    }

    double AAPsiCompMove(double tuning, int nrep)    {

        double nacc = 0;
        double ntot = 0;

        for (int rep=0; rep<nrep; rep++)    {
            double delta = -GetLogPrior();
            double m = tuning * (Random::Uniform() - 0.5);
            double e = exp(m);
            psi /= e;
            for (int i=0; i<Naarr; i++) {
                aadist[i] *= e;
            }
            delta += GetLogPrior();
            delta += (Naarr-1)*m;

            int accepted = (log(Random::Uniform()) < delta);
            if (accepted) {
                nacc++;
            } else {
                psi *= e;
                for (int i=0; i<Naarr; i++) {
                    aadist[i] /= e;
                }
            }
            ntot++;
        }

        Update();
        return nacc / ntot;
    }

    double GetMeanAAEntropy() const {
        double m1 = 0;
        double m0 = 0;
        for (int i=0; i<Naa; i++)   {
            for (int j=0; j<Gcat; j++)  {
                m0 += obsweight[i] * Gweight[j];
                m1 += obsweight[i] * Gweight[j] * Random::GetEntropy(selacprofiles->GetVal(i,j));
            }
        }
        return m1 / m0;
    }

    double GetMeanAADist() const {
        double m1 = 0;
        for (int i=0; i<Naarr; i++) {
            m1 += aadist[i];
        }
        m1 /= Naarr;
        return m1;
    }

    double GetVarAADist() const {
        double m1 = 0;
        double m2 = 0;
        for (int i=0; i<Naarr; i++) {
            m1 += aadist[i];
            m2 += aadist[i]*aadist[i];
        }
        m1 /= Naarr;
        m2 /= Naarr;
        m2 -= m1*m1;
        return m2;
    }

    void TraceHeader(ostream &os) const override {
        os << "#logprob\tGinvshape\tpsi\tmeanaa\tvaraa\taaent\tomega\tvaromega\n";
    }

    void Trace(ostream& os) const override {
        os << GetLogProb() << '\t';
        os << Ginvshape << '\t';
        os << psi << '\t';
        os << GetMeanAADist() << '\t';
        os << GetVarAADist() << '\t';
        os << GetMeanAAEntropy() << '\t';
        os << predmeanomega << '\t';
        os << predrelvaromega << '\n';
    }

    // miscell.

    void PostPred(string name) override {}

    void Monitor(ostream &os) const override {
        os << "aa fitness\t";
        for (size_t i=0; i<aadist_acc.size(); i++)  {
            os << aadist_acc[i] / aadist_tot[i] << '\t';
        }
        os << '\n';

        os << "G\t";
        for (size_t i=0; i<G_acc.size(); i++)  {
            os << G_acc[i] / G_tot[i] << '\t';
        }
        os << '\n';

        os << "psi\t";
        for (size_t i=0; i<psi_acc.size(); i++)  {
            os << psi_acc[i] / psi_tot[i] << '\t';
        }
        os << '\n';

        os << "G/psi\t";
        for (size_t i=0; i<gpsi_acc.size(); i++)  {
            os << gpsi_acc[i] / gpsi_tot[i] << '\t';
        }
        os << '\n';
    }

    void MonitorCompAAFreqs(ostream& os) const  {
        double tot = 0;
        for (int i=0; i<Naarr; i++)   {
            tot += aadist[i];
        }
        tot /= Naarr;
        for (int a=0; a<Naa; a++)   {
            for (int b=0; b<Naa; b++)   {
                os << AminoAcids[a] << AminoAcids[b] << '\t' << predcompaafreq[a][b] << '\t' << obscompaafreq[a][b] << '\t' << log(predcompaafreq[a][b]) << '\t' << log(obscompaafreq[a][b]) << '\t';
               if (a == b)  {
                  os << 0 << '\n';
               }
               else {
                  os << aadist[rrindex(a,b)]/tot << '\n';
               }
            }
        }
    }

    void MonitorAADist(ostream& os) {
        double tot = 0;
        for (int i=0; i<Naarr; i++)   {
            tot += aadist[i];
        }
        tot /= Naarr;
        int index = 0;
        for (int a=0; a<Naa; a++)   {
            for (int b=a+1; b<Naa; b++) {
                os << AminoAcids[a] << AminoAcids[b] << '\t' << aadist[index]/tot << '\t' << granthaadist[index] << '\n';
                index++;
            }
        }
    }

    void FromStream(istream &is) override {
        is >> Ginvshape;
        is >> psi;
        for (int i=0; i<Naarr; i++)   {
            is >> aadist[i];
        }
    }

    void ToStream(ostream &os) const override {
        os << Ginvshape;
        double tot = 0;
        for (int i=0; i<Naarr; i++)   {
            tot += aadist[i];
        }
        tot /= Naarr;
        os << '\t' << psi * tot;
        for (int i=0; i<Naarr; i++)   {
            os << '\t' << aadist[i] / tot;
        }
        os << '\n';
    }


    void Run()   {
        ofstream tos((name + ".trace").c_str());
        ofstream aos((name + ".chain").c_str());
        TraceHeader(tos);
        while(1)    {
            Move();
            Trace(tos);
            tos.flush();
            ToStream(aos);
            aos.flush();
            ofstream mos((name + ".monitor").c_str());
            Monitor(mos);
            ofstream fos((name + ".aafreqs").c_str());
            MonitorCompAAFreqs(fos);
            ofstream aos((name + ".aadist").c_str());
            MonitorAADist(aos);
        }
    }
};

int main(int argc, char* argv[])    {

    if (argc == 1) {
        cerr << "fastselac -kappa <val> -gc <val> -omega <mean> <relvar> -u <val> -epsilon <val> <name>\n";
        exit(1);
    }

    double gamma = 0.5;
    double kappa = 3;
    double obsmeanomega = 0.3;
    double obsrelvaromega = 2;
    double epsilon = 1.0;
    double u = 0;
    int Gcat = 4;
    string aafitnessfile = "";
    string name = "";

    int i = 1;
    while (i < argc) {
        string s = argv[i];

        if (s == "-omega") {
            i++;
            obsmeanomega = atof(argv[i]);
            i++;
            obsrelvaromega = atof(argv[i]);
        }
        else if (s == "-kappa") {
            i++;
            kappa = atof(argv[i]);
        }
        else if (s == "-gc")    {
            i++;
            gamma = atof(argv[i]);
        }
        else if (s == "-u") {
            i++;
            u = atof(argv[i]);
        }
        else if (s == "-eps")   {
            i++;
            epsilon = atof(argv[i]);
        }
        else if (s == "-gcat")   {
            i++;
            Gcat = atoi(argv[i]);
        }
        else if (s == "-aa")    {
            i++;
            aafitnessfile = argv[i];
        }
        else    {
            if (i != argc-1)    {
                cerr << "error in command\n";
                exit(1);
            }
            name = argv[i];
        }
        i++;
    }

    cerr << "create model\n";
    FastEmpiricalSelACModel model(aafitnessfile, Gcat, kappa, gamma, u, obsmeanomega, obsrelvaromega, epsilon, name);
    cerr << "run model\n";
    model.Run();
}



