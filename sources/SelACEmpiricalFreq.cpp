
#include "CodonSequenceAlignment.hpp"
#include "Random.hpp"
#include <fstream>

int main(int argc, char* argv[])    {


    vector<vector<double>> empfreq(Naa, vector<double>(Naa, 0));
    vector<double> weight(Naa, 0);

    ifstream is(argv[1]);
    ofstream os(argv[2]);

    int Ngene;
    is >> Ngene;
    for (int g=0; g<Ngene; g++) {
        string datafile;
        is >> datafile;
        cerr << datafile << '\n';
        auto data = new FileSequenceAlignment(datafile);
        if (data->GetNsite() % 3) {
            cerr << "error : not a correctly formatted codon-alignment: " << datafile << '\n';
            exit(1);
        }
        auto codondata = new CodonSequenceAlignment(data, true);
        auto count = codondata->GetSiteAAEmpiricalCounts();
        for (size_t i=0; i<count.size(); i++)   {
            int max = 0;
            int jmax = 0;
            for (size_t j=0; j<Naa; j++)    {
                if (count[i][j] > max)  {
                    max = count[i][j];
                    jmax = j;
                }
            }
            if (max)    {
                weight[jmax]++;
                for (size_t j=0; j<Naa; j++)    {
                    empfreq[jmax][j] += count[i][j];
                }
            }
        }
    }

    for (int i=0; i<Naa; i++)   {
        double tot = 0;
        for (int j=0; j<Naa; j++)   {
            tot += empfreq[i][j];
        }
        for (int j=0; j<Naa; j++)   {
            empfreq[i][j] /= tot;
        }
    }

    double tot = 0;
    for (int i=0; i<Naa; i++)   {
        tot += weight[i];
    }
    for (int i=0; i<Naa; i++)   {
        weight[i] /= tot;
    }

    os << Naa << '\n';
    for (int i=0; i<Naa; i++)   {
        os << weight[i];
        for (int j=0; j<Naa; j++)   {
            os << '\t' << empfreq[i][j];
        }
        os << '\n';
    }
}

