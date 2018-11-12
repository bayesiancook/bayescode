#include <cmath>
#include <fstream>
#include "AAMutSelDSBDPOmegaModel.hpp"
#include "components/ChainDriver.hpp"
#include "components/ChainReader.hpp"
#include "components/stats_posterior.hpp"
#include "tclap/CmdLine.h"

using namespace std;
using namespace TCLAP;

class ReadAAMutSelDSBDPOmegaArgParse {
    CmdLine &cmd;

  public:
    ReadAAMutSelDSBDPOmegaArgParse(CmdLine &cmd) : cmd(cmd) {}

    ValueArg<int> every{
        "e", "every", "Number of iterations between two traces", false, 1, "int", cmd};
    ValueArg<int> until{"u", "until", "Maximum number of (saved) iterations (-1 means unlimited)",
        false, -1, "int", cmd};
    ValueArg<int> burnin{"b", "burnin", "Number of iterations for burnin", false, 0, "int", cmd};
    SwitchArg ppred{"p", "ppred",
        "For each point of the chain (after burn-in), produces a data replicate simulated "
        "from the posterior predictive distribution",
        cmd};
    SwitchArg ss{
        "s", "ss", "Computes the mean posterior site-specific state equilibrium frequencies", cmd};
    UnlabeledValueArg<std::string> chain_name{
        "chain_name", "Chain name (output file prefix)", true, "chain", "string", cmd};
};

int main(int argc, char *argv[]) {
    CmdLine cmd{"AAMutSelDSBDPOmega", ' ', "0.1"};
    ReadAAMutSelDSBDPOmegaArgParse read_args(cmd);
    cmd.parse(argc, argv);

    int burnin = read_args.burnin.getValue();
    int every = read_args.every.getValue();
    int until = read_args.until.getValue();
    // pveber until is set to -1, meaning to unlimited as described
    int size = (until - burnin) / every;
    std::string chain_name = read_args.chain_name.getValue();

    std::ifstream is{chain_name + ".param"};
    ChainDriver::fake_read(is);  // We're not interested in the ChainDriver of the param file
    AAMutSelDSBDPOmegaModel model(is);
    ChainReader cr{model, chain_name + ".chain"};

    cr.skip(burnin);
    cerr << size << " points to read\n";

    if (read_args.ppred.getValue()) {
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
        stats_posterior<AAMutSelDSBDPOmegaModel>(model, cr, every, size);
    }
}