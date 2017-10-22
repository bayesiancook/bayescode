#include <cmath>
#include <fstream>
#include "Chain.hpp"
#include "DiscGammaSiteOmegaModel.hpp"
using namespace std;

class DiscGammaSiteOmegaChain : public Chain  {

  private:
    // Chain parameters
    string modeltype, datafile, treefile;
    int Ncat;

  public:
    DiscGammaSiteOmegaModel* GetModel() {
        return static_cast<DiscGammaSiteOmegaModel*>(model);
    }

    string GetModelType() override { return modeltype; }

    DiscGammaSiteOmegaChain(string indatafile, string intreefile, int inNcat, int inevery, int inuntil, string inname, int force) : modeltype("DISCGAMMASITEOMEGA"), datafile(indatafile), treefile(intreefile), Ncat(inNcat) {
        every = inevery;
        until = inuntil;
        name = inname;
        New(force);
    }

    DiscGammaSiteOmegaChain(string filename) {
        name = filename;
        Open();
        Save();
    }

    void New(int force) override {
        model = new DiscGammaSiteOmegaModel(datafile,treefile,Ncat);
        GetModel()->Allocate();
        GetModel()->Unfold();
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
        is >> Ncat;
        int tmp;
        is >> tmp;
        if (tmp) {
            cerr << "-- Error when reading model\n";
            exit(1);
        }
        is >> every >> until >> size;

        if (modeltype == "DISCGAMMASITEOMEGA") {
            model = new DiscGammaSiteOmegaModel(datafile,treefile,Ncat);
        } else {
            cerr << "-- Error when opening file " << name
                 << " : does not recognise model type : " << modeltype << '\n';
            exit(1);
        }
        GetModel()->Allocate();
        model->FromStream(is);
        model->Update();
        GetModel()->Unfold();
        cerr << size << " points saved, current ln prob = " << GetModel()->GetLogProb() << "\n";
        model->Trace(cerr);
    }

    void Save() override {
        ofstream param_os((name + ".param").c_str());
        param_os << GetModelType() << '\n';
        param_os << datafile << '\t' << treefile << '\n';
        param_os << Ncat << '\n';
        param_os << 0 << '\n';
        param_os << every << '\t' << until << '\t' << size << '\n';
        model->ToStream(param_os);
    }
};

int main(int argc, char* argv[])	{

    // starting a chain from existing files
    if (argc == 2 && argv[1][0] != '-') {
        string name = argv[1];
        DiscGammaSiteOmegaChain* chain = new DiscGammaSiteOmegaChain(name);
        cerr << "chain " << name << " started\n";
        chain->Start();
        cerr << "chain " << name << " stopped\n";
        cerr << chain->GetSize() << " points saved, current ln prob = " << chain->GetModel()->GetLogProb() << "\n";
        chain->GetModel()->Trace(cerr);
    }

    // new chain
    else    {
        string datafile = "";
        string treefile = "";
        int Ncat = 4;
        string name = "";
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
                else if (s == "-dgam") {
                    i++;
                    Ncat = atoi(argv[i]);
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
            cerr << "discom -d <alignment> -t <tree> -dgam <ncat> <chainname> \n";
            cerr << '\n';
            exit(1);
        }

        DiscGammaSiteOmegaChain* chain = new DiscGammaSiteOmegaChain(datafile,treefile,Ncat,every,until,name,force);
        cerr << "chain " << name << " started\n";
        chain->Start();
        cerr << "chain " << name << " stopped\n";
        cerr << chain->GetSize() << "-- Points saved, current ln prob = " << chain->GetModel()->GetLogProb() << "\n";
        chain->GetModel()->Trace(cerr);
    }
}


