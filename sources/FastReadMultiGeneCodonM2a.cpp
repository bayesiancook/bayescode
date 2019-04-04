
#include <cmath>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include <algorithm>

using namespace std;

int PostAnalysis(string newpath, string chain_name, int burnin, int with_sites)   {

    cerr << chain_name << '\n';
    cerr << "get parameters\n";
    // open param file, get ali file or list
    ifstream is((chain_name + ".param").c_str());
    string model_name, data_path, data_name, tree_file;
    is >> model_name >> data_path >> data_name >> tree_file;

    if (newpath != "None")  {
        data_path = newpath;
    }

    cerr << "get gene list\n";
    int Ngene;
    ifstream genelist_is((chain_name + ".genelist").c_str());
    genelist_is >> Ngene;
    vector<string> gene_list(Ngene,"");
    map<string,int> gene_nsite;
    for (int gene=0; gene<Ngene; gene++) {
        genelist_is >> gene_list[gene] >> gene_nsite[gene_list[gene]];
    }

    // open ali file or list and get original gene list and gene nsite's
    ifstream data_is((data_path + data_name).c_str());
    string tmp;
    data_is >> tmp;
    vector<string> original_gene_list(Ngene,"");
    if (tmp == "ALI")   {
        int ngene;
        data_is >> ngene;
        if (ngene != Ngene) {
            cerr << "error: non matching number of genes\n";
            cerr << ngene << '\t' << Ngene << '\n';
            exit(1);
        }
        for (int gene=0; gene<Ngene; gene++) {
            data_is >> original_gene_list[gene];
            int ntaxa, nsite;
            data_is >> ntaxa >> nsite;
            for (int j=0; j<ntaxa; j++) {
                data_is >> tmp >> tmp;
            }
        }
    }
    else    {
        cerr << data_path << '\n';
        cerr << data_name << '\n';
        cerr << tmp << '\n';
        cerr << "implement for lists\n";
        exit(1);
    }

    cerr << "process hyper params\n";
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
    for (int k=0; k<nvar; k++)  {
        mean[k] /= size;
    }

    const string hyperparam[] = {"purom_mean", "purom_invconc", "dposom_mean", "dposom_invshape", "purw_mean", "purw_invconc", "posw_mean", "posw_invconc", "pi"};
    ofstream hyper_os((chain_name + ".posthyper").c_str());
    for (int k=0; k<nvar; k++)  {
        hyper_os << hyperparam[k] << '\t' << mean[k] << '\t' << min[k] << '\t' << max[k] << '\n';
    }

    cerr << "process gene pps\n";
    // open .genepp file and get post mean gene pps
    ifstream posw_is((chain_name + ".posw").c_str());
    for (int i=0; i<burnin; i++)    {
        getline(posw_is,line);
    }

    int check_size = 0;
    map<string, double> posw;
    map<string, double> pp;
    for (int gene=0; gene<Ngene; gene++)    {
        posw[gene_list[gene]] = 0;
        pp[gene_list[gene]] = 0;
    }
    while (getline(posw_is,line))  {
        istringstream is(line);
        check_size++;
        for (int gene=0; gene<Ngene; gene++)    {
            double tmp;
            is >> tmp;
            posw[gene_list[gene]] += tmp;
            if (tmp)    {
                pp[gene_list[gene]] ++;
            }
        }
    }
    if (check_size != size) {
        cerr << "error: file size does not match (posw)\n";
        exit(1);
    }
    for (int gene=0; gene<Ngene; gene++)    {
        posw[gene_list[gene]] /= size;
        pp[gene_list[gene]] /= size;
    }

    ifstream posom_is((chain_name + ".posom").c_str());
    for (int i=0; i<burnin; i++)    {
        getline(posom_is,line);
    }

    check_size = 0;
    map<string, double> posom;
    map<string, double> min_posom;
    map<string, double> max_posom;
    map<string, vector<double>> gene_posom_sample;
    for (int gene=0; gene<Ngene; gene++)    {
        posom[gene_list[gene]] = 0;
        gene_posom_sample[gene_list[gene]].assign(size,0);
    }
    while (getline(posom_is,line))  {
        istringstream is(line);
        for (int gene=0; gene<Ngene; gene++)    {
            double tmp;
            is >> tmp;
            posom[gene_list[gene]] += tmp;
            gene_posom_sample[gene_list[gene]][check_size] = tmp;
        }
        check_size++;
    }
    if (check_size != size) {
        cerr << "error: file size does not match (posw)\n";
        exit(1);
    }
    for (int gene=0; gene<Ngene; gene++)    {
        posom[gene_list[gene]] /= size;
        sort(gene_posom_sample[gene_list[gene]].begin(), gene_posom_sample[gene_list[gene]].end());
        min_posom[gene_list[gene]] = gene_posom_sample[gene_list[gene]][(int) (alpha * size)];
        max_posom[gene_list[gene]] = gene_posom_sample[gene_list[gene]][(int) ((1-alpha) * size)];
    }


    map<string, vector<double>> site_pp;
    if (with_sites) {
        cerr << "process site pps\n";
        // open .sitepp file and get site pps

        ifstream sitepp_is((chain_name + ".sitepp").c_str());
        for (int i=0; i<burnin; i++)    {
            getline(sitepp_is,line);
        }

        for (int gene=0; gene<Ngene; gene++)    {
            site_pp[gene_list[gene]].assign(gene_nsite[gene_list[gene]], 0);
        }

        check_size = 0;
        while (getline(sitepp_is,line))  {
            istringstream is(line);
            for (int gene=0; gene<Ngene; gene++)    {
                string name;
                is >> name;
                if (name != gene_list[gene])    {
                    cerr << "error when reading site pps: non matching gene name\n";
                    cerr << name << " instead of " << gene_list[gene] << '\n';
                    exit(1);
                }
                int nsite = gene_nsite[name];
                for (int i=0; i<nsite; i++) {
                    double tmp;
                    is >> tmp;
                    site_pp[name][i] += tmp;
                }
            }
            check_size ++;
        }
        if (check_size != size) {
            cerr << "error: file size does not match (sitepp)\n";
            exit(1);
        }
        for (int gene=0; gene<Ngene; gene++)    {
            string name = gene_list[gene];
            for (int i=0; i<gene_nsite[name]; i++)  {
                site_pp[name][i] /= size;
            }
        }
    }

    cerr << "write output\n";
    // write gene and site post pps
    ofstream post_os((chain_name + ".postanalysis").c_str());
    for (int gene=0; gene<Ngene; gene++)    {
        string name = original_gene_list[gene];
        post_os << name << '\t' << pp[name] << '\t' << posw[name] << '\t' << posom[name] << '\t' << min_posom[name] << '\t' << max_posom[name];
        if (with_sites) {
            for (int i=0; i<gene_nsite[name]; i++) {
                post_os << '\t' << (int) (100 * site_pp[name][i]);
            }
        }
        post_os << '\n';
    }
    return 1;
}

int main(int argc, char *argv[]) {

    string newpath = "None";
    int burnin = 0;
    int with_sites = 0;
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
            } else if (s == "-s")   {
                with_sites = 0;
            } else if (s == "+s")   {
                with_sites = 1;
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
        cerr << "readglobom [-x <burnin> +s/-s] <chainname> \n";
        cerr << '\n';
        exit(1);
    }
    PostAnalysis(newpath, name, burnin, with_sites);
}
