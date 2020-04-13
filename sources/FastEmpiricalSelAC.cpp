#include "Random.hpp"
#include "grantham.hpp"
#include "LG.hpp"
#include "MeanSubMatrix.hpp"
#include "HKYSubMatrix.hpp"
#include "EmpiricalSelACProfileBidimArray.hpp"
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
    // vector<int> aamult;

    // amino acid selection
    int Gcat;
    double Ginvshape;
    vector<double> G;
    vector<double> Gweight;
    double psi;

    vector<double> obsweight;
    vector<vector<double>> obscompaafreq;
    vector<vector<double>> predcompaafreq;

    vector<vector<double>> aafitness;

    // putting things together
    EmpiricalSelACProfileBidimArray* selacprofiles;

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

    vector<double> aafit_acc;
    vector<double> aafit_tot;
    vector<double> G_acc;
    vector<double> G_tot;
    vector<double> psi_acc;
    vector<double> psi_tot;

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
        aafit_acc.assign(1,0);
        aafit_tot.assign(1,0.1);
        G_acc.assign(3,0);
        G_tot.assign(3,0.1);
        psi_acc.assign(3,0);
        psi_tot.assign(3,0.1);

        nucmatrix = new HKYSubMatrix(kappa,gamma);

        Ginvshape = 2.0;
        G.assign(Gcat, 1.0/Gcat);
        UpdateG();
        psi = 3.5;
        statespace = new CodonStateSpace(Universal);
        /*
        aamult.assign(Naa, 0);
        for (int i=0; i<statespace->GetNstate(); i++)   {
            int aa = statespace->Translation(i);
            aamult[aa]++;
        }
        for (int aa=0; aa<Naa; aa++)    {
            cerr << AminoAcids[aa] << '\t' << aamult[aa] << '\n';
        }
        */

        aafitness.assign(Naa, vector<double>(Naa,1.0));
        for (int i=0; i<Naa; i++)   {
            for (int j=0; j<Naa; j++)   {
                aafitness[i][j] = obscompaafreq[i][j];
            }
        }

        Gweight.assign(Gcat, 1.0/Gcat);

        selacprofiles = new EmpiricalSelACProfileBidimArray(aafitness,G,psi);

        codonmatrices = new AAMutSelPolyMutCodonMatrixBidimArray(*selacprofiles, *statespace, *nucmatrix, u, 1.0);

        Update();
    }

    // updates and logprob 

    void UpdateG()  {
        Random::DiscGamma(G,1.0/Ginvshape);
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

    void Update(int a)  {
        selacprofiles->UpdateRow(a);
        codonmatrices->CorruptRow(a);
        ComputePredictedDNDS(a);
        ComputeSummaryStatistics();
    }

    double GetMSD() const {
        double diff = 0;
        double tmp = 100*(obsmeanomega - predmeanomega)/obsmeanomega;
        diff += tmp*tmp;
        tmp = 100*(obsrelvaromega - predrelvaromega)/obsrelvaromega;
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
        for (int i=0; i<Naa; i++)   {
            for (int j=0; j<Naa; j++)   {
                tot -= aafitness[i][j];
            }
        }
        return tot;
    }

    double GetLogPrior() const {
        double tot = 0;
        tot -= log(psi);
        tot -= log(Ginvshape);
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
            aafit_acc[0] += MoveAAFitness(200, 1.0);
            aafit_tot[0] ++;
            MoveAllAAFitness(10, 1.0);
            MoveGinvshape();
            // MovePsi();
        }
        return 1.0;
    }

    double MoveAAFitness(int nrep, double tuning) {
        double nacc = 0;
        double ntot = 0;
        for (int i=0; i<nrep; i++)  {
            int a = (int) (Naa * Random::Uniform());
            int b = (int) (Naa * Random::Uniform());
            double delta = -GetLogProb();
            double m = tuning * (Random::Uniform() - 0.5);
            double e = exp(m);
            aafitness[a][b] *= e;
            Update(a);
            delta += GetLogProb();
            delta += m;
            int acc = (log(Random::Uniform()) < delta);
            if (acc)    {
                nacc++;
            }
            else    {
                aafitness[a][b] /= e;
                Update(a);
            }
            ntot++;
        }
        return nacc / ntot;
    }

    double MoveAllAAFitness(int nrep, double tuning)    {
        double nacc = 0;
        double ntot = 0;
        for (int i=0; i<nrep; i++)  {
            int a = (int) (Naa * Random::Uniform());
            double delta = -GetLogPrior();
            double m = tuning * (Random::Uniform() - 0.5);
            double e = exp(m);
            for (int b=0; b<Naa; b++)   {
                aafitness[a][b] *= e;
            }
            Update(a);
            delta += GetLogPrior();
            delta += m;
            int acc = (log(Random::Uniform()) < delta);
            if (acc)    {
                nacc++;
            }
            else    {
                for (int b=0; b<Naa; b++)   {
                    aafitness[a][b] /= e;
                }
                Update(a);
            }
            ntot++;
        }
        return nacc / ntot;
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

    void TraceHeader(ostream &os) const override {
        os << "#logprob\tGinvshape\tpsi\taaent\tomega\tvaromega\n";
    }

    void Trace(ostream& os) const override {
        os << GetLogProb() << '\t';
        os << Ginvshape << '\t';
        os << psi << '\t';
        os << GetMeanAAEntropy() << '\t';
        os << predmeanomega << '\t';
        os << predrelvaromega << '\n';
    }

    // miscell.

    void PostPred(string name) override {}

    void Monitor(ostream &os) const override {
        os << "aa fitness\t";
        for (size_t i=0; i<aafit_acc.size(); i++)  {
            os << aafit_acc[i] / aafit_tot[i] << '\t';
        }
        os << '\n';

        os << "G\t";
        for (size_t i=0; i<G_acc.size(); i++)  {
            os << G_acc[i] / G_tot[i] << '\t';
        }
        os << '\n';

        os << "psi\t";
        for (size_t i=0; i<G_acc.size(); i++)  {
            os << psi_acc[i] / psi_tot[i] << '\t';
        }
        os << '\n';
    }

    void MonitorCompAAFreqs(ostream& os) const  {
        for (int a=0; a<Naa; a++)   {
            double tot = 0;
            for (int b=0; b<Naa; b++)   {
                tot += aafitness[a][b];
            }
            for (int b=0; b<Naa; b++)   {
                os << AminoAcids[a] << AminoAcids[b] << '\t' << predcompaafreq[a][b] << '\t' << obscompaafreq[a][b] << '\t' << log(predcompaafreq[a][b]) << '\t' << log(obscompaafreq[a][b]) << '\t' << aafitness[a][b] / tot * Naa << '\t' << log(aafitness[a][b] / tot * Naa) << '\n';
            }
        }
    }

    void FromStream(istream &is) override {
        is >> Ginvshape;
        is >> psi;
        for (int i=0; i<Naa; i++)   {
            for (int j=0; j<Naa; j++)   {
                is >> aafitness[i][j];
            }
        }
    }

    void ToStream(ostream &os) const override {
        os << Ginvshape;
        os << '\t' << psi;
        for (int i=0; i<Naa; i++)   {
            double tot = 0;
            for (int j=0; j<Naa; j++)   {
                tot += aafitness[i][j];
            }
            for (int j=0; j<Naa; j++)   {
                os << '\t' << aafitness[i][j] / tot * Naa;
            }
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
    double epsilon = 0.001;
    double u = 0.01;
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
        else if (s == "-cat")   {
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



