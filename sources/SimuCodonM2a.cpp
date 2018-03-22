#include <cmath>
#include <fstream>
#include "CodonM2aModel.hpp"
using namespace std;


int main(int argc, char* argv[])	{

    string datafile = "";
    string treefile = "";
    string paramfile = "";
    string name = "";

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
            else if (s == "-p") {
                i++;
                paramfile = argv[i];
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
        cerr << "simucodonm2a -d <alignment> -t <tree> -p <paramfile> <name> \n";
        cerr << '\n';
        exit(1);
    }

    double pi = 0.1;
    cerr << "new model\n";
    CodonM2aModel* model = new CodonM2aModel(datafile,treefile,pi);
    model->Allocate();
    cerr << "set lengths\n";
    model->SetLengthsFromTree();
    ifstream param_is(paramfile.c_str());
    model->FromStreamCodeML(param_is);
    cerr << "post pred\n";
    model->PostPred(name);
}



