#include <cmath>
#include <fstream>
#include "Chain.hpp"
#include "AAMutSelSparseM9Model.hpp"
using namespace std;

/**
 * \brief Chain object for running an MCMC under AAMutSelSparseM9Model 
 */

class AAMutSelSparseM9Chain : public Chain  {

  private:
    // Chain parameters
    string modeltype, datafile, treefile;
    int fixhyper;
    int fixfitness;
    // -1: estimated
    double epsilon;
    double pi;
    double pospi;
    double posmeanhypermean, posmeanhyperinvshape;
    double posinvshapehypermean, posinvshapehyperinvshape;
    double poswhypermean, poswhyperinvconc;

  public:
    AAMutSelSparseM9Model* GetModel() {
        return static_cast<AAMutSelSparseM9Model*>(model);
    }

    string GetModelType() override { return modeltype; }

    AAMutSelSparseM9Chain(string indatafile, string intreefile, 
        int infixhyper, int infixfitness, double inepsilon, double inpi, 
        double inpospi, 
        double inposmeanhypermean, double inposmeanhyperinvshape,
        double inposinvshapehypermean, double inposinvshapehyperinvshape,
        double inposwhypermean, double inposwhyperinvconc,
        int inevery, int inuntil, string inname, int force) : 
        
        modeltype("AAMUTSELSPARSEM9"), datafile(indatafile), treefile(intreefile), 
        fixhyper(infixhyper), fixfitness(infixfitness), epsilon(inepsilon), pi(inpi),
        pospi(inpospi),
        posmeanhypermean(inposmeanhypermean),
        posmeanhyperinvshape(inposmeanhyperinvshape),
        posinvshapehypermean(inposinvshapehypermean),
        posinvshapehyperinvshape(inposinvshapehyperinvshape),
        poswhypermean(inposwhypermean),
        poswhyperinvconc(inposwhyperinvconc)    {

        every = inevery;
        until = inuntil;
        name = inname;
        New(force);
    }

    AAMutSelSparseM9Chain(string filename) {
        name = filename;
        Open();
        Save();
    }

    void New(int force) override {
        cerr << "new model\n";
        model = new AAMutSelSparseM9Model(datafile,treefile,fixhyper,fixfitness,epsilon,pi,pospi);

        GetModel()->SetMixtureHyperParameters(
                pospi,
                posmeanhypermean, posmeanhyperinvshape,
                posinvshapehypermean, posinvshapehyperinvshape,
                poswhypermean, poswhyperinvconc);

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
        is >> fixhyper;
        is >> fixfitness;
        is >> epsilon;
        is >> pi;
        is >> pospi;
        is >> posmeanhypermean >> posmeanhyperinvshape;
        is >> posinvshapehypermean >> posinvshapehyperinvshape;
        is >> poswhypermean >> poswhyperinvconc;
        int tmp;
        is >> tmp;
        if (tmp) {
            cerr << "-- Error when reading model\n";
            exit(1);
        }
        is >> every >> until >> size;

        if (modeltype == "AAMUTSELSPARSEM9") {
            model = new AAMutSelSparseM9Model(datafile,treefile,fixhyper,fixfitness,epsilon,pi,pospi);
        } else {
            cerr << "-- Error when opening file " << name
                 << " : does not recognise model type : " << modeltype << '\n';
            exit(1);
        }
        GetModel()->SetChainSize(GetSize());

        GetModel()->SetMixtureHyperParameters(
                pospi,
                posmeanhypermean, posmeanhyperinvshape,
                posinvshapehypermean, posinvshapehyperinvshape,
                poswhypermean, poswhyperinvconc);

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
        param_os << fixhyper << '\n';
        param_os << fixfitness << '\n';
        param_os << epsilon << '\n';
        param_os << pi << '\n';
        param_os << pospi << '\n';
        param_os << posmeanhypermean << '\t' << posmeanhyperinvshape << '\n';
        param_os << posinvshapehypermean << '\t' << posinvshapehyperinvshape << '\n';
        param_os << poswhypermean << '\t' << poswhyperinvconc << '\n';
        param_os << 0 << '\n';
        param_os << every << '\t' << until << '\t' << size << '\n';
        model->ToStream(param_os);
    }
};

int main(int argc, char* argv[])	{

    string name = "";
    AAMutSelSparseM9Chain* chain = 0;

    // starting a chain from existing files
    if (argc == 2 && argv[1][0] != '-') {
        name = argv[1];
        chain = new AAMutSelSparseM9Chain(name);
    }

    // new chain
    else    {
        string datafile = "";
        string treefile = "";
        int fixhyper = 3;
        int fixfitness = 0;
        double epsilon = -1;
        double pi = -1;

        double pospi = 0.1;
        double poswhypermean = 0.5;
        double poswhyperinvconc = 0.1;

        double posmeanhypermean = 1.0;
        double posmeanhyperinvshape = 1.0;
        double posinvshapehypermean = 1.0;
        double posinvshapehyperinvshape = 1.0;

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
                else if (s == "-pospi") {
                    i++;
                    pospi = atof(argv[i]);
                }
                else if (s == "-dposom") {
                    i++;
                    posmeanhypermean = atof(argv[i]);
                    i++;
                    posmeanhyperinvshape = atof(argv[i]);
                    i++;
                    posinvshapehypermean = atof(argv[i]);
                    i++;
                    posinvshapehyperinvshape = atof(argv[i]);
                } 
                else if (s == "-posw") {
                    i++;
                    poswhypermean = atof(argv[i]);
                    i++;
                    poswhyperinvconc = atof(argv[i]);
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
            cerr << "aamutselsparsem9 -d <alignment> -t <tree> -ncat <ncat> <chainname> \n";
            cerr << '\n';
            exit(1);
        }

        chain = new AAMutSelSparseM9Chain(datafile,treefile,
                fixhyper,fixfitness,epsilon,pi,
                pospi,
                posmeanhypermean, posmeanhyperinvshape, 
                posinvshapehypermean, posinvshapehyperinvshape, 
                poswhypermean, poswhyperinvconc, 
                every,until,name,force);
    }

    cerr << "chain " << name << " started\n";
    chain->Start();
    cerr << "chain " << name << " stopped\n";
    cerr << chain->GetSize() << "-- Points saved, current ln prob = " << chain->GetModel()->GetLogProb() << "\n";
    chain->GetModel()->Trace(cerr);
}


