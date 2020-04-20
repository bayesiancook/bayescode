#include <cmath>
#include <fstream>
#include "Chain.hpp"
#include "AAMutSelSparseOmegaModel.hpp"
using namespace std;

/**
 * \brief Chain object for running an MCMC under AAMutSelSparseOmegaModel 
 */

class AAMutSelSparseOmegaChain : public Chain  {

  private:
    // Chain parameters
    string modeltype, datafile, treefile;
    int omegamode, omegaprior;
    double dposompi, dposomhypermean, dposomhyperinvshape;
    double maxdposom;
    int fixhyper;
    int fixfitness;
    // -1: estimated
    double epsilon;
    double pi;

  public:
    AAMutSelSparseOmegaModel* GetModel() {
        return static_cast<AAMutSelSparseOmegaModel*>(model);
    }

    string GetModelType() override { return modeltype; }

    AAMutSelSparseOmegaChain(string indatafile, string intreefile, 
        int inomegamode, int inomegaprior,
        double indposompi, double indposomhypermean,
        double indposomhyperinvshape, double inmaxdposom,
        int infixhyper, int infixfitness, double inepsilon, double inpi, 
        int inevery, int inuntil, string inname, int force) : 
        
        modeltype("AAMUTSELSparseOMEGA"), datafile(indatafile), treefile(intreefile), 
        omegamode(inomegamode),
        omegaprior(inomegaprior),
        dposompi(indposompi),
        dposomhypermean(indposomhypermean),
        dposomhyperinvshape(indposomhyperinvshape),
        maxdposom(inmaxdposom),
        fixhyper(infixhyper), fixfitness(infixfitness), epsilon(inepsilon), pi(inpi)    {

        every = inevery;
        until = inuntil;
        name = inname;
        New(force);
    }

    AAMutSelSparseOmegaChain(string filename) {
        name = filename;
        Open();
        Save();
    }

    void New(int force) override {
        cerr << "new model\n";
        model = new AAMutSelSparseOmegaModel(datafile,treefile,omegamode,omegaprior,fixhyper,fixfitness,epsilon,pi);
        if (omegaprior != 0) {
            GetModel()->SetMaxDPosOm(maxdposom);

            GetModel()->SetDPosOmHyperParameters(dposompi, dposomhypermean, dposomhyperinvshape);
        }

        cerr << "allocate\n";
        GetModel()->Allocate();
        cerr << "update\n";
        GetModel()->Update();
        cerr << "-- Reset" << endl;
        Reset(force);
        cerr << "-- initial ln prob = " << GetModel()->GetLogProb() << "\n";
        model->Trace(cerr);
    }

    void Open() override {
        ifstream is((name + ".param").c_str());
        if (!is) {
            cerr << "-- Error : cannot find file : " << name << ".param\n";
            exit(1);
        }
        is >> modeltype;
        is >> datafile >> treefile;
        is >> omegamode >> omegaprior >> dposompi >> dposomhypermean >> dposomhyperinvshape >> maxdposom;
        is >> fixhyper;
        is >> fixfitness;
        is >> epsilon;
        is >> pi;
        int tmp;
        is >> tmp;
        if (tmp) {
            cerr << "-- Error when reading model\n";
            exit(1);
        }
        is >> every >> until >> size;

        if (modeltype == "AAMUTSELSparseOMEGA") {
            model = new AAMutSelSparseOmegaModel(datafile,treefile,omegamode,omegaprior,fixhyper,fixfitness,epsilon,pi);
        } else {
            cerr << "-- Error when opening file " << name
                 << " : does not recognise model type : " << modeltype << '\n';
            exit(1);
        }
        if (omegaprior != 0) {
            GetModel()->SetMaxDPosOm(maxdposom);
            GetModel()->SetDPosOmHyperParameters(dposompi, dposomhypermean,
                                                 dposomhyperinvshape);
        }
        GetModel()->SetChainSize(GetSize());
        GetModel()->Allocate();
        model->FromStream(is);
        GetModel()->Update();
        cerr << size << " points saved, current ln prob = " << GetModel()->GetLogProb() << "\n";
        model->Trace(cerr);
    }

    void Save() override {
        ofstream param_os((name + ".param").c_str());
        param_os << GetModelType() << '\n';
        param_os << datafile << '\t' << treefile << '\n';
        param_os << omegamode << '\t' << omegaprior;
        param_os << '\t' << dposompi << '\t' << dposomhypermean;
        param_os << '\t' << dposomhyperinvshape << '\t' << maxdposom << '\n';
        param_os << fixhyper << '\n';
        param_os << fixfitness << '\n';
        param_os << epsilon << '\n';
        param_os << pi << '\n';
        param_os << 0 << '\n';
        param_os << every << '\t' << until << '\t' << size << '\n';
        model->ToStream(param_os);
    }
};

int main(int argc, char* argv[])	{

    string name = "";
    AAMutSelSparseOmegaChain* chain = 0;

    // starting a chain from existing files
    if (argc == 2 && argv[1][0] != '-') {
        name = argv[1];
        chain = new AAMutSelSparseOmegaChain(name);
    }

    // new chain
    else    {
        string datafile = "";
        string treefile = "";
        int omegamode = 3;
        int omegaprior = 0;
        double dposompi = 0.1;
        double dposomhypermean = 1.0;
        double dposomhyperinvshape = 0.5;
        double maxdposom = 0;
        int fixhyper = 3;
        int fixfitness = 0;
        double epsilon = -1;
        double pi = -1;
        name = "";
        int force = 1;
        int every = 1;
        int until = -1;

        try	{

            if (argc == 1)	{
                throw(0);
            }

            int i = 1;
            while (i < argc)	{
                string s = argv[i];

                if (s == "-d")	{
                    i++;
                    datafile = argv[i];
                }
                else if ((s == "-t") || (s == "-T"))	{
                    i++;
                    treefile = argv[i];
                }
                else if (s == "-f")	{
                    force = 1;
                } 
                else if (s == "-maxdposom")   {
                    i++;
                    maxdposom = atof(argv[i]);
                } 
                else if (s == "-fixomega") {
                    omegamode = 3;
                } 
                else if (s == "-freeomega") {
                    omegamode = 1;
                } 
                else if (s == "-postomega") {
                    omegamode = -1;
                }
                else if (s == "-gamomega") {
                    omegaprior = 0;
                } 
                else if ((s == "-mixomega") || (s == "-gammixomega")) {
                    omegaprior = 1;
                    i++;
                    dposompi = atof(argv[i]);
                    i++;
                    dposomhypermean = atof(argv[i]);
                    i++;
                    dposomhyperinvshape = atof(argv[i]);
                } 
                else if (s == "-loggammixomega") {
                    omegaprior = 2;
                    i++;
                    dposompi = atof(argv[i]);
                    i++;
                    dposomhypermean = atof(argv[i]);
                    i++;
                    dposomhyperinvshape = atof(argv[i]);
                } 
                else if (s == "-cauchymixomega") {
                    omegaprior = 3;
                    i++;
                    dposompi = atof(argv[i]);
                    dposomhypermean = 0;
                    i++;
                    dposomhyperinvshape = atof(argv[i]);
                }
                else if (s == "-fixhyper")  {
                    fixhyper = 3;
                }
                else if (s == "-freehyper") {
                    fixhyper = 0;
                }
                else if ((s == "-eps") || (s == "-epsilon"))   {
                    i++;
                    string tmp = argv[i];
                    if (tmp == "free")  {
                        epsilon = -1;
                    }
                    else    {
                        epsilon = atof(argv[i]);
                    }
                }
                else if (s == "-pi")    {
                    i++;
                    pi = atof(argv[i]);
                }
                else if (s == "-fixaa") {
                    fixfitness = 1;
                }
                else if ( (s == "-x") || (s == "-extract") )	{
                    i++;
                    if (i == argc) throw(0);
                    every = atoi(argv[i]);
                    i++;
                    if (i == argc) throw(0);
                    until = atoi(argv[i]);
                }
                else	{
                    if (i != (argc -1))	{
                        throw(0);
                    }
                    name = argv[i];
                }
                i++;
            }
            if ((datafile == "") || (treefile == "") || (name == ""))	{
                throw(0);
            }
        }
        catch(...)	{
            cerr << "aamutseldp -d <alignment> -t <tree> -ncat <ncat> <chainname> \n";
            cerr << '\n';
            exit(1);
        }

        chain = new AAMutSelSparseOmegaChain(datafile,treefile,
                omegamode, omegaprior, dposompi,
                dposomhypermean, dposomhyperinvshape, maxdposom, 
                fixhyper,fixfitness,epsilon,pi,every,until,name,force);
    }

    cerr << "chain " << name << " started\n";
    chain->Start();
    cerr << "chain " << name << " stopped\n";
    cerr << chain->GetSize() << "-- Points saved, current ln prob = " << chain->GetModel()->GetLogProb() << "\n";
    chain->GetModel()->Trace(cerr);
}


