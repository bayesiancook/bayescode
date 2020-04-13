
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

class FastSelACModel : public ProbModel {

    public:

    string name;

    double epsilon;

    double kappa;
    double gamma;
    HKYSubMatrix *nucmatrix;

    double u;

    CodonStateSpace* statespace;

    // amino acid selection
    int Gcat;
    double Ginvshape;
    vector<double> G;
    vector<double> Gweight;
    double psi;

    vector<double> aadist;
    vector<double> aaweight;

    // putting things together
    SelACProfileBidimArray* selacprofiles;
    // to do with double and triple mutants
    AAMutSelPolyMutCodonMatrixBidimArray* codonmatrices;
    MeanCodonSubMatrixFromBidimArray* meancodonmatrix;
    AAProjectedCodonSubMatrix* aamatrix;

    // 20 amino acid frequencies
    vector<double> obsaafreq;
    vector<double> predaafreq;

    // array of 75 fixation biases
    vector<double> obsaarelrate;
    vector<double> predaarelrate;

    // mean dN/dS across sites
    double obsmeanomega;
    double predmeanomega;

    // variance in dN/dS across sites
    double obsrelvaromega;
    double predrelvaromega;

    int Naa;
    int Naarr;

    vector<double> aadist_acc;
    vector<double> aadist_tot;
    vector<double> aaweight_acc;
    vector<double> aaweight_tot;
    vector<double> G_acc;
    vector<double> G_tot;
    vector<double> psi_acc;
    vector<double> psi_tot;

    // allocation

    FastSelACModel(int inGcat, double inkappa, double ingamma, double inu, double inobsmeanomega, double inobsrelvaromega, double inepsilon, string inname) {

        name = inname;

        Gcat = inGcat;
        kappa = inkappa;
        gamma = ingamma;
        u = inu;
        obsmeanomega = inobsmeanomega;
        obsrelvaromega = inobsrelvaromega;
        epsilon = inepsilon;

        Naa = 20;
        Naarr = Naa * (Naa-1) / 2;

        obsaafreq.assign(Naa,0);
        predaafreq.assign(Naa,0);
        for (int i=0; i<Naa; i++)   {
            obsaafreq[i] = LG_Stat[i];
        }

        obsaarelrate.assign(Naarr,0);
        predaarelrate.assign(Naarr,0);
        for (int i=0; i<Naarr; i++)   {
            obsaarelrate[i] = LG_RR[i];
        }

        // monitoring
        aadist_acc.assign(3,0);
        aadist_tot.assign(3,0.1);
        aaweight_acc.assign(3,0);
        aaweight_tot.assign(3,0.1);
        G_acc.assign(3,0);
        G_tot.assign(3,0.1);
        psi_acc.assign(3,0);
        psi_tot.assign(3,0.1);

        // allocate
        nucmatrix = new HKYSubMatrix(kappa,gamma);
        Ginvshape = 2.0;
        G.assign(Gcat, 1.0/Gcat);
        UpdateG();
        psi = 1.0;
        statespace = new CodonStateSpace(Universal);

        aadist.assign(Naarr, 1.0);
        /*
        vector<double> center(Naarr, 1.0/Naarr);
        double conc = Naarr;
        Random::DirichletSample(aadist,center,conc);
        */
        UpdateGrantham(grantham_wcom_selac, grantham_wpol, grantham_wvol);

        aaweight.assign(Naa, 1.0/Naa);
        for (int i=0; i<Naa; i++)   {
            aaweight[i] = LG_Stat[i];
        }

        Gweight.assign(Gcat, 1.0/Gcat);

        selacprofiles = new SelACProfileBidimArray(aadist,G,psi);
        codonmatrices = new AAMutSelPolyMutCodonMatrixBidimArray(*selacprofiles, *statespace, *nucmatrix, u, 1.0);
        meancodonmatrix = new MeanCodonSubMatrixFromBidimArray(*statespace, *codonmatrices, aaweight, Gweight); 
        aamatrix = new AAProjectedCodonSubMatrix(meancodonmatrix);
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

    void UpdateMatrices() {
        /*
        nucmatrix->SetKappa(kappa);
        nucmatrix->SetGC(gamma);
        */
        codonmatrices->Corrupt();
        meancodonmatrix->ComputeFullArray();
        aamatrix->ComputeFullArray();
    }

    void GetPredictedDNDS(double& mean, double& relvar) const  {
        double m0 = 0;
        double m1 = 0;
        double m2 = 0;
        for (int i=0; i<Naa; i++) {
            for (int j=0; j<Gcat; j++)  {
                double dnds = codonmatrices->GetVal(i,j).GetPredictedDNDS();
                m0 += aaweight[i] * Gweight[j];
                m1 += aaweight[i] * Gweight[j] * dnds;
                m2 += aaweight[i] * Gweight[j] * dnds * dnds;
            }
        }
        m1 /= m0;
        m2 /= m0;
        m2 -= m1*m1;
        mean = m1;
        relvar = m2 / m1 / m1;
    }

    void ComputeSummaryStatistics() {
        GetPredictedDNDS(predmeanomega, predrelvaromega);
        aamatrix->GetStationaryArray(predaafreq);
        aamatrix->GetRelRateArray(predaarelrate);
    }

    void Update() override {
        UpdateG();
        selacprofiles->SetPsi(psi);
        selacprofiles->Update();
        UpdateMatrices();
        ComputeSummaryStatistics();
    }

    void Update(int a, int b)  {

        selacprofiles->UpdateRow(a);
        selacprofiles->UpdateRow(b);
        meancodonmatrix->DenormalizeByStat();
        meancodonmatrix->AddRow(a,-1);
        meancodonmatrix->AddRow(b,-1);
        codonmatrices->CorruptRow(a);
        codonmatrices->CorruptRow(b);
        meancodonmatrix->AddRow(a,+1);
        meancodonmatrix->AddRow(b,+1);
        meancodonmatrix->NormalizeByStat();

        aamatrix->ComputeFullArray();
        ComputeSummaryStatistics();
    }

    double GetMSD() const {
        double diff = 0;
        double tmp = 100*(obsmeanomega - predmeanomega)/obsmeanomega;
        diff += tmp*tmp;
        tmp = 100*(obsrelvaromega - predrelvaromega)/obsrelvaromega;
        diff += tmp*tmp;
        for (int i=0; i<Naa; i++)   {
            tmp = 20 * (obsaafreq[i] - predaafreq[i]);
            diff += tmp*tmp;
        }
        double obstot = 0;
        double predtot = 0;
        for (int i=0; i<Naarr; i++)   {
            obstot += obsaarelrate[i];
            predtot += predaarelrate[i];
            tmp = obsaarelrate[i] - predaarelrate[i];
            diff += tmp*tmp;
        }
        return diff;
    }

    double GetLogPrior() const {
        double tot = 0;
        for (int i=0; i<Naarr; i++) {
            tot -= aadist[i];
        }
        tot -= log(psi);
        tot -= log(Ginvshape);
        return tot;
    }

    double GetLogProb() const override {
        double diff = GetMSD();
        return GetLogPrior() - diff * diff / 2 / epsilon / epsilon;
    }

    // moving

    double Move() override {
        aadist_acc[0] += FastMoveAADist(200,1.0);
        aadist_tot[0] ++;
        // MoveAADist();
        // MoveAAWeight();
        // MoveGinvshape();
        MovePsi();
        AAPsiCompMove(1,10);
        AAPsiCompMove(0.1,10);
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

    double MoveAADist() {
        int nrep = 50;
        for (int i=0; i<nrep; i++) {
            // int j = (int) (Naarr * Random::Uniform());
            int a = (int) (Naa * Random::Uniform());
            int b = (int) ((Naa-1) * Random::Uniform());
            if (b >= a) b++;
            int j = rrindex(a,b);
            aadist_acc[0] += ScalingMove(aadist[j], 1.00, 1, &FastSelACModel::GetLogProb, &FastSelACModel::Update, this);
            aadist_tot[0] ++;
        }
        /*
        aadist_acc[0] += ProfileMove(aadist, 1.00, 1, 10, &FastSelACModel::GetLogProb, &FastSelACModel::Update, this);
        aadist_tot[0] ++;
        aadist_acc[1] += ProfileMove(aadist, 0.3, 3, 10, &FastSelACModel::GetLogProb, &FastSelACModel::Update, this);
        aadist_tot[1] ++;
        aadist_acc[2] += ProfileMove(aadist, 0.1, 3, 10, &FastSelACModel::GetLogProb, &FastSelACModel::Update, this);
        aadist_tot[2] ++;
        */
        return 1.0;
    }

    double MoveAAWeight() {
        aaweight_acc[0] += ProfileMove(aaweight, 1.00, 1, 10, &FastSelACModel::GetLogProb, &FastSelACModel::Update, this);
        aaweight_tot[0] ++;
        aaweight_acc[1] += ProfileMove(aaweight, 0.3, 3, 10, &FastSelACModel::GetLogProb, &FastSelACModel::Update, this);
        aaweight_tot[1] ++;
        aaweight_acc[2] += ProfileMove(aaweight, 0.1, 3, 10, &FastSelACModel::GetLogProb, &FastSelACModel::Update, this);
        aaweight_tot[2] ++;
        return 1.0;
    }

    double MoveGinvshape()  {
        G_acc[0] += ScalingMove(Ginvshape, 0.3, 3, &FastSelACModel::GetLogProb, &FastSelACModel::Update, this);
        G_tot[0] ++;
        G_acc[1] += ScalingMove(Ginvshape, 0.03, 3, &FastSelACModel::GetLogProb, &FastSelACModel::Update, this);
        G_tot[1] ++;
        G_acc[2] += ScalingMove(Ginvshape, 0.003, 3, &FastSelACModel::GetLogProb, &FastSelACModel::Update, this);
        G_tot[2] ++;
        return 1.0;
    }

    double MovePsi()  {
        psi_acc[0] += ScalingMove(psi, 0.3, 3, &FastSelACModel::GetLogProb, &FastSelACModel::Update, this);
        psi_tot[0] ++;
        psi_acc[1] += ScalingMove(psi, 0.03, 3, &FastSelACModel::GetLogProb, &FastSelACModel::Update, this);
        psi_tot[1] ++;
        psi_acc[2] += ScalingMove(psi, 0.003, 3, &FastSelACModel::GetLogProb, &FastSelACModel::Update, this);
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
                m0 += aaweight[i] * Gweight[j];
                m1 += aaweight[i] * Gweight[j] * Random::GetEntropy(selacprofiles->GetVal(i,j));
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
        os << "#logprob\tGinvshape\tpsi\tmeanaa\tvaraa\taaweight\taaent\tomega\tvaromega\n";
    }

    void Trace(ostream& os) const override {
        os << GetLogProb() << '\t';
        os << Ginvshape << '\t';
        os << psi << '\t';
        os << GetMeanAADist() << '\t';
        os << GetVarAADist() << '\t';
        os << Random::GetEntropy(aaweight) << '\t';
        os << GetMeanAAEntropy() << '\t';
        os << predmeanomega << '\t';
        os << predrelvaromega << '\n';
    }

    // miscell.

    void PostPred(string name) override {}

    void Monitor(ostream &os) const override {
        os << "aadist\t";
        for (size_t i=0; i<aadist_acc.size(); i++)  {
            os << aadist_acc[i] / aadist_tot[i] << '\t';
        }
        os << '\n';

        os << "aaweight\t";
        for (size_t i=0; i<aaweight_acc.size(); i++)  {
            os << aaweight_acc[i] / aaweight_tot[i] << '\t';
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

    void MonitorRelRates(ostream& os) const  {
        int i = 0;
        for (int a=0; a<Naa; a++)   {
            for (int b=a+1; b<Naa; b++) {
                os << AminoAcids[a] << AminoAcids[b] << '\t' << predaarelrate[i] << '\t' << obsaarelrate[i] << '\n';
                i++;
            }
        }
    }

    void MonitorAAFreqs(ostream& os) const  {
        for (int a=0; a<Naa; a++)   {
            os << AminoAcids[a] << '\t' << predaafreq[a] << '\t' << obsaafreq[a] << '\n';
        }
    }

    void MonitorAADist(ostream &os) const {
        for (int i=0; i<Naarr; i++)   {
            os << aadist[i] << '\n';
        }
    }

    void MonitorAAWeights(ostream &os) const {
        for (int i=0; i<Naa; i++)   {
            os << aaweight[i] << '\n';
        }
        for (int i=0; i<Naarr; i++)   {
            os << '\t' << aadist[i];
        }
        os << '\n';
    }

    void FromStream(istream &is) override {
        is >> Ginvshape;
        for (int i=0; i<Naa; i++)   {
            is >> aaweight[i];
        }
        for (int i=0; i<Naarr; i++)   {
            is >> aadist[i];
        }
    }

    void ToStream(ostream &os) const override {
        os << Ginvshape;
        for (int i=0; i<Naa; i++)   {
            os << '\t' << aaweight[i];
        }
        for (int i=0; i<Naarr; i++)   {
            os << '\t' << aadist[i];
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
            ofstream ros((name + ".rr").c_str());
            MonitorRelRates(ros);
            ofstream fos((name + ".aafreqs").c_str());
            MonitorAAFreqs(fos);
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
    FastSelACModel model(Gcat, kappa, gamma, u, obsmeanomega, obsrelvaromega, epsilon, name);
    cerr << "run model\n";
    model.Run();
}



