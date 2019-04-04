
#include <cmath>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <map>

using namespace std;

int PostAnalysis(string newpath, string chain_name, int burnin, int every, int until)   {

    // open param file, get ali file or list
    ifstream is((chain_name + ".param").c_str());
    string model_name, data_path, data_name, tree_file;
    is >> model_name >> data_path >> data_name >> tree_file;

    if (newpath != "None")  {
        data_path = newpath;
    }

    int Ngene;
    ifstream genelist_is((chain_name + ".genelist").c_str());
    is >> Ngene;
    vector<string> gene_list(Ngene,"");
    for (int gene=0; gene<Ngene; gene++) {
        is >> gene_list[gene];
    }

    // open ali file or list and get original gene list and gene nsite's
    ifstream data_is((data_path + data_name).c_str());
    string tmp;
    is >> tmp;
    vector<string> original_gene_list(Ngene,"");
    map<string,int> gene_nsite;
    if (tmp == "ALI")   {
        int ngene;
        is >> ngene;
        if (ngene != Ngene) {
            cerr << "error: non matching number of genes\n";
            exit(1);
        }
        for (int gene=0; gene<Ngene; gene++) {
            is >> original_gene_list[gene];
            int ntaxa, nsite;
            is >> ntaxa >> nsite;
            gene_nsite[original_gene_list[gene]] = nsite;
            for (int j=0; j<ntaxa; j++) {
                is >> tmp >> tmp;
            }
        }
    }
    else    {
        cerr << "implement for lists\n";
        exit(1);
    }

    // open .chain file and get post mean hyperparameters
    ifstream chain_is((chain_name + ".chain").c_str());
    string line;
    for (int i=0; i<burnin; i++)    {
        getline(chain_is,line);
    }

    int nvar = 9;
    int size = 0;
    vector<vector<double>> sample(9,vector<double>(0,0));
    while (getline(chain_is,line))  {
        istringstream is(line);
        for (int k=0; k<nvar; k++) {
            double s;
            is >> s;
            sample[k].push_back(s);
        }
        size++;
    }
    vector<double> mean(nvar,0);
    vector<double> min(nvar,0);
    vector<double> max(nvar,0);
    for (int i=0; i<size; i++)  {
        for (int k=0; k<nvar; k++)  {
            mean[k] += sample[k][i];
        }
    }
    double alpha = 0.025;
    for (int k=0; k<nvar; k++)  {
        sort(sample[k].begin(), sample[k].end());
        min[k] = sample[k][(int) (alpha * size)];
        max[k] = sample[k][(int) ((1-alpha) * size)];
    }

    const string hyperparam[] = {"purom_mean", "purom_invconc", "dposom_mean", "dposom_invshape", "purw_mean", "purw_invconc", "posw_mean", "posw_invconc", "pi"};
    ofstream hyper_os((chain_name + ".posthyper").c_str());
    for (int k=0; k<nvar; k++)  {
        hyper_os << hyperparam[k] << '\t' << mean[k] << '\t' << min[k] << '\t' << max[k] << '\n';
    }

    // open .genepp file and get post mean gene pps
    ifstream posw_is((chain_name + ".posw").c_str());
    for (int i=0; i<burnin; i++)    {
        getline(posw_is,line);
    }

    int check_size = 0;
    vector<double> posw(Ngene,0);
    while (getline(posw_is,line))  {
        istringstream is(line);
        check_size++;
        for (int gene=0; gene<Ngene; gene++)    {
            double tmp;
            is >> tmp;
            posw[gene] += tmp;
        }
    }
    if (check_size != size) {
        cerr << "error: file size does not match (posw)\n";
        exit(1);
    }
    for (int gene=0; gene<Ngene; gene++)    {
        posw[gene] /= size;
    }

    ifstream posom_is((chain_name + ".posom").c_str());
    for (int i=0; i<burnin; i++)    {
        getline(posom_is,line);
    }

    check_size = 0;
    vector<double> posom(Ngene,0);
    while (getline(posom_is,line))  {
        istringstream is(line);
        check_size++;
        for (int gene=0; gene<Ngene; gene++)    {
            double tmp;
            is >> tmp;
            posom[gene] += tmp;
        }
    }
    if (check_size != size) {
        cerr << "error: file size does not match (posw)\n";
        exit(1);
    }
    for (int gene=0; gene<Ngene; gene++)    {
        posom[gene] /= size;
    }


    // open .sitepp file and get site pps

    // write gene and site post pps
    ofstream post_os((chain_name + ".postanalysis").c_str());
    for (int gene=0; gene<Ngene; gene++)    {
    }
    return 1;
}

int main(int argc, char *argv[]) {

    string newpath = "None";
    int burnin = 0;
    int every = 1;
    int until = -1;
    string name;

    try {
        if (argc == 1) {
            throw(0);
        }

        int i = 1;
        while (i < argc) {
            string s = argv[i];
            if ((s == "-x") || (s == "-extract")) {
                i++;
                burnin = atoi(argv[i]);
            } else if (s == "-p")   {
                i++;
                newpath = argv[i];
            } else {
                if (i != (argc - 1)) {
                    throw(0);
                }
                name = argv[i];
            }
            i++;
        }
        if (name == "") {
            throw(0);
        }
    } catch (...) {
        cerr << "readglobom [-x <burnin> <every> <until>] <chainname> \n";
        cerr << '\n';
        exit(1);
    }
    PostAnalysis(newpath, name, burnin, every, until);
}
