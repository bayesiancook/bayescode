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

    SwitchArg ss{"s", "ss",
        "Computes the mean posterior site-specific amino-acid equilibrium frequencies"
        "(amino-acid fitness profiles). "
        "Results are written in {chain_name}.siteprofiles by default (optionally use the --output "
        "argument "
        " to specify a different output path).",
        cmd};
    ValueArg<string> omega_pp{"", "omega_threshold",
        "Threshold to compute the mean posterior probability that ω⁎ "
        "(or ω if option `flatfitness` is used in `mutselomega`) is greater than a given value "
        "(1.0 to test for adaptation). "
        "Results are written in {chain_name}.omegappgt{omega_pp}.tsv by default (optionally use "
        "the --output argument to specify a different output path).",
        false, "", "string", cmd};
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
    } else if (read_args.trace.getValue()) {
        string file_name = read_args.OutputFile(".trace.tsv");
        recompute_trace<AAMutSelDM5Model>(model, cr, file_name, every, size);
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

        string file_name = read_args.OutputFile(".siteprofiles");
        ofstream os(file_name);
        os << "site\tA\tC\tD\tE\tF\tG\tH\tI\tK\tL\tM\tN\tP\tQ\tR\tS\tT\tV\tW\tY\n";
        for (int i = 0; i < model.GetNsite(); i++) {
            os << i + 1;
            for (auto &aa : sitestat[i]) {
                aa /= size;
                os << '\t' << aa;
            }
            os << '\n';
        }
        cerr << "mean site-specific profiles in " << file_name << "\n";
        cerr << '\n';
    } else if (!read_args.omega_pp.getValue().empty()) {
        double omega_pp = stod(read_args.omega_pp.getValue());
        std::vector<double> omegappgto(model.GetNsite(), 0);
        std::vector<double> omega(model.GetNsite(), 0);

        for (int step = 0; step < size; step++) {
            cerr << '.';
            cr.skip(every);
            for (int site = 0; site < model.GetNsite(); site++) {
                omega[site] += model.GetSiteOmega(site);
                if (model.GetSiteOmega(site) > omega_pp) { omegappgto[site]++; }
            }
        }
        cerr << '\n';

        string filename =
            read_args.OutputFile(".omegappgt" + read_args.omega_pp.getValue() + ".tsv");
        ofstream os(filename);
        for (int i = 0; i < model.GetNsite(); i++) {
            os << i + 1 << '\t' << omegappgto[i] / size << '\t' << omega[i] / size << '\n';
        }
        cerr << "Posterior prob of omega greater than " << read_args.omega_pp.getValue() << " in "
             << filename << "\n";
        cerr << '\n';
    } else {
        stats_posterior<AAMutSelDM5Model>(model, cr, every, size);
    }
}