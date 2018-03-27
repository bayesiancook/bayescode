#include <cmath>
#include <fstream>
#include "Chain.hpp"
#include "AAMutSelDSBDPOmegaNeModel.hpp"
using namespace std;

/**
 * \brief Chain object for running an MCMC under AAMutSelDSBDPOmegaNeModel
 */

class AAMutSelDSBDPOmegaNeChain : public Chain  {

  private:
    // Chain parameters
    string modeltype, datafile, treefile;
    int Ncat, baseNcat, Ncond;

  public:
    AAMutSelDSBDPOmegaNeModel* GetModel() {
        return static_cast<AAMutSelDSBDPOmegaNeModel*>(model);
    }

    string GetModelType() override { return modeltype; }

    AAMutSelDSBDPOmegaNeChain(string indatafile, string intreefile, int inNcat, int inbaseNcat, int inNcond, int inevery, int inuntil, string inname, int force) : modeltype("AAMUTSELDSBDPOMEGANE"), datafile(indatafile), treefile(intreefile), Ncat(inNcat), baseNcat(inbaseNcat), Ncond(inNcond) {
        every = inevery;
        until = inuntil;
        name = inname;
        New(force);
    }

    AAMutSelDSBDPOmegaNeChain(string filename) {
        name = filename;
        Open();
        Save();
    }

    void New(int force) override {
        cerr << "new model\n";
        model = new AAMutSelDSBDPOmegaNeModel(datafile,treefile,Ncat,baseNcat,Ncond);
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
        is >> Ncat >> baseNcat >> Ncond;
        int tmp;
        is >> tmp;
        if (tmp) {
            cerr << "-- Error when reading model\n";
            exit(1);
        }
        is >> every >> until >> size;

        if (modeltype == "AAMUTSELDSBDPOMEGANE") {
            model = new AAMutSelDSBDPOmegaNeModel(datafile,treefile,Ncat,baseNcat,Ncond);
        } else {
            cerr << "-- Error when opening file " << name
                 << " : does not recognise model type : " << modeltype << '\n';
            exit(1);
        }
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
        param_os << Ncat << '\t' << baseNcat << '\t' << Ncond << '\n';
        param_os << 0 << '\n';
        param_os << every << '\t' << until << '\t' << size << '\n';
        model->ToStream(param_os);
    }
};

int main(int argc, char* argv[])	{

    string name = "";
    AAMutSelDSBDPOmegaNeChain* chain = 0;

    // starting a chain from existing files
    if (argc == 2 && argv[1][0] != '-') {
        name = argv[1];
        chain = new AAMutSelDSBDPOmegaNeChain(name);
    }

    // new chain
    else    {
        string datafile = "";
        string treefile = "";
        int Ncat = 100;
        int baseNcat = 1;
        int Ncond = 2;
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
                else if (s == "-ncat") {
                    i++;
                    Ncat = atoi(argv[i]);
                }
                else if (s == "-basencat") {
                    i++;
                    baseNcat = atoi(argv[i]);
                }
                else if (s == "-ncond") {
                    i++;
                    Ncond = atoi(argv[i]);
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
            cerr << "aamutselddpne -d <alignment> -t <tree> -ncat <ncat> -ncond <ncond> <chainname> \n";
            cerr << '\n';
            exit(1);
        }

        chain = new AAMutSelDSBDPOmegaNeChain(datafile,treefile,Ncat,baseNcat,Ncond,every,until,name,force);
    }

    cerr << "chain " << name << " started\n";
    chain->Start();
    cerr << "chain " << name << " stopped\n";
    cerr << chain->GetSize() << "-- Points saved, current ln prob = " << chain->GetModel()->GetLogProb() << "\n";
    chain->GetModel()->Trace(cerr);
}
