#include <cmath>
#include <fstream>
#include "AAMutSelDM5Model.hpp"
#include "components/ChainDriver.hpp"
#include "components/ChainReader.hpp"
#include "components/ReadArgParse.hpp"
#include "components/stats_posterior.hpp"
#include "tclap/CmdLine.h"

using namespace std;
using namespace TCLAP;

class ReadAAMutSelDSBDPOmegaArgParse : public ReadArgParse {
  public:
    explicit ReadAAMutSelDSBDPOmegaArgParse(CmdLine &cmd) : ReadArgParse(cmd) {}

    SwitchArg ss{
        "s", "ss", "Computes the mean posterior site-specific state equilibrium frequencies", cmd};
    ValueArg<double> omega_pp{
        "", "omega_threshold", "the threshold for omega", false, 1.0, "double", cmd};
};

int main(int argc, char *argv[]) {
    CmdLine cmd{"AAMutSelDM5", ' ', "0.1"};
    ReadAAMutSelDSBDPOmegaArgParse read_args(cmd);
    cmd.parse(argc, argv);

    std::string chain_name = read_args.GetChainName();
    int burnin = read_args.GetBurnIn();
    int every = read_args.GetEvery();
    int size = read_args.GetSize();

    std::ifstream is{chain_name + ".param"};
    ChainDriver::fake_read(is);  // We're not interested in the ChainDriver of the param file
    AAMutSelDM5Model model(is);
    ChainReader cr{model, chain_name + ".chain"};

    cr.skip(burnin);
    cerr << size << " points to read\n";

    if (read_args.GetPpred()) {
        for (int i = 0; i < size; i++) {
            cerr << '.';
            cr.skip(every);
            model.PostPred("ppred_" + chain_name + "_" + std::to_string(i) + ".ali");
        }
        cerr << '\n';
    } else if (read_args.ss.getValue()) {
        std::vector<std::vector<double>> sitestat(model.GetNsite(), {0});

        for (int step = 0; step < size; step++) {
            cerr << '.';
            cr.skip(every);
            for (int i = 0; i < model.GetNsite(); i++) {
                std::vector<double> const &profile = model.GetProfile(i);
                if (sitestat[i].size() != profile.size()) {
                    sitestat[i].resize(profile.size(), 0);
                };
                for (unsigned k{0}; k < profile.size(); k++) { sitestat[i][k] += profile[k]; }
            }
        }
        cerr << '\n';

        ofstream os((chain_name + ".siteprofiles").c_str());
        os << model.GetNsite() << '\n';
        for (int i = 0; i < model.GetNsite(); i++) {
            os << i + 1;
            for (auto &aa : sitestat[i]) {
                aa /= size;
                os << '\t' << aa;
            }
            os << '\n';
        }
        cerr << "mean site-specific profiles in " << chain_name << ".siteprofiles\n";
        cerr << '\n';
    } else {
        std::vector<double> omegappgto(model.GetNsite(), 0);
        std::vector<double> omega(model.GetNsite(), 0);

        for (int step = 0; step < size; step++) {
            cerr << '.';
            cr.skip(every);
            for (int site = 0; site < model.GetNsite(); site++) {
                omega[site] += model.GetSiteOmega(site);
                if (model.GetSiteOmega(site) > read_args.omega_pp.getValue()) {
                    omegappgto[site]++;
                }
            }
        }
        cerr << '\n';

        string filename{chain_name + ".omegappgt" + to_string(read_args.omega_pp.getValue())};
        ofstream os(filename.c_str());
        for (int i = 0; i < model.GetNsite(); i++) {
            os << i + 1 << '\t' << omegappgto[i] / size << '\t' << omega[i] / size << '\n';
        }
        cerr << "Posterior prob of omega greater than " << read_args.omega_pp.getValue() << " in "
             << filename << "\n";
        cerr << '\n';
    }
}