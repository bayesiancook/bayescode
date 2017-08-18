#include <cmath>
#include <fstream>
#include "CodonM2aModel.hpp"
using namespace std;

int main(int argc, char* argv[])	{

    string datafile = "";
    string treefile = "";
    double pi = 0.1;
    string paramname = "";
    string hyperparamname = "";
    string name = "";
    // int force = 1;

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
            /*
            else if (s == "-f")	{
                force = 1;
            }
            */
            else if (s == "-pi")    {
                i++;
                pi = atof(argv[i]);
            }
            else if (s == "-fixparam")  {
                i++;
                paramname = argv[i];
            }
            else if (s == "-fixhyper")  {
                i++;
                hyperparamname = argv[i];
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
        cerr << "codonm2a -d <alignment> -t <tree> [-fixparam <paramfile> -fixhyper <hyperparamfile>] <chainname> \n";
        cerr << '\n';
        exit(1);
    }

	CodonM2aModel* model = new CodonM2aModel(datafile,treefile,pi);
    model->Allocate();
    /*
    if (paramname != "")    {
        ifstream is(paramname.c_str());
        model->GetGlobalParametersFromFile(is);
        model->FixGlobalParameters();
        model->GetHyperParametersFromFile(is);
    }
    if (hyperparamname != "")    {
        ifstream is(hyperparamname.c_str());
        model->GetHyperParametersFromFile(is);
    }
    */
    model->Unfold();
	ofstream os((name + ".trace").c_str());
	ofstream pos((name + ".sitepp").c_str());
	model->TraceHeader(os);
	os.flush();
	while(1)	{
		model->Move();
		model->Trace(os);
        model->TracePostProb(pos);
		os.flush();
        pos.flush();
	}
}

