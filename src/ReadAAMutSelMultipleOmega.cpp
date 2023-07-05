#include <cmath>
#include <fstream>
#include "AAMutSelMultipleOmegaModel.hpp"
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

    SwitchArg nuc{"n", "nuc",
        "Mean posterior 4x4 nucleotide matrix."
        "Results are written in {chain_name}.nucmatrix.tsv.",
        cmd};
    ValueArg<string> chain_omega{"", "chain_omega",
        "A second chain ran with the option --freeomega and --flatfitness to obtain the classical "
        "ω-based codon model (Muse & Gaut). "
        "These two chains allow to compute posterior of ω, ω₀, ωᴬ=ω-ω₀ and p(ωᴬ>0) for each site "
        "and at the gene level. "
        "Results are written in {chain_name}.omegaA.tsv.",
        false, "", "string", cmd};
    SwitchArg omega{"", "omega",
        "Compute posterior credible interval for ω for each site and at the gene level. "
        "Can be combined with the option `confidence_interval` to change the default value (0.025 "
        "at each side of the distribution). "
        "Results are written in {chain_name}.ci{confidence_interval}.tsv.",
        cmd};
    SwitchArg omega_knot{"", "omega_0",
        "Compute posterior credible interval for ω₀ predicted at the mutation-selection "
        "equilibrium from the fitness profiles, for each site and at the gene level. "
        "Can be combined with the option `confidence_interval` to change the default value (0.025 "
        "at each side of the distribution). "
        "Results are written in {chain_name}.ci{confidence_interval}.tsv.",
        cmd};
    ValueArg<string> confidence_interval{"c", "confidence_interval",
        "Boundary for posterior credible interval of ω and ω₀ (per site and at the gene level). "
        "Default value is 0.025 at each side, meaning computing the 1-2*0.025=95% CI.",
        false, "0.025", "string", cmd};
    SwitchArg ss{"s", "ss",
        "Computes the mean posterior site-specific amino-acid equilibrium frequencies"
        "(amino-acid fitness profiles). "
        "Results are written in {chain_name}.siteprofiles.",
        cmd};
    TCLAP::ValueArg<string> profiles{"", "profiles",
        "Change the profiles filename if desired, "
        "otherwise given by {chain_name}.siteprofiles as default.", false,
        "", "string", cmd};
    ValueArg<string> omega_pp{"", "omega_threshold",
        "Threshold to compute the mean posterior probability that ω⁎ "
        "(or ω if option `flatfitness` is used in `mutselomega`) is greater than a given value. "
        "Results are written in {chain_name}.omegappgt{omega_pp}.tsv.",
        false, "1.0", "string", cmd};

    string GetProfilesName() {
        if (profiles.getValue().empty()) {
            return GetChainName() + ".siteprofiles";
        } else {
            return profiles.getValue();
        }
    }
};

int main(int argc, char *argv[]) {
    CmdLine cmd{"AAMutSelMultipleOmega", ' ', "0.1"};
    ReadAAMutSelDSBDPOmegaArgParse read_args(cmd);
    cmd.parse(argc, argv);

    string chain_name = read_args.GetChainName();
    int burnin = read_args.GetBurnIn();
    int every = read_args.GetEvery();
    int size = read_args.GetSize();

    ifstream is{chain_name + ".param"};
    ChainDriver::fake_read(is);  // We're not interested in the ChainDriver of the param file
    AAMutSelMultipleOmegaModel model(is);
    ChainReader cr{model, chain_name + ".chain"};

    cr.skip(burnin);
    cerr << size << " points to read\n";

    if (read_args.GetPpred()) {
        for (int i = 0; i < size; i++) {
            cerr << '.';
            cr.skip(every);
            model.PostPred("ppred_" + chain_name + "_" + to_string(i) + ".ali");
        }
        cerr << '\n';
    } else if (read_args.ss.getValue()) {
        vector<vector<double>> sitestat(model.GetNsite(), {0});

        for (int step = 0; step < size; step++) {
            cerr << '.';
            cr.skip(every);
            for (int i = 0; i < model.GetNsite(); i++) {
                vector<double> const &profile = model.GetProfile(i);
                if (sitestat[i].size() != profile.size()) {
                    sitestat[i].resize(profile.size(), 0);
                };
                for (unsigned k{0}; k < profile.size(); k++) { sitestat[i][k] += profile[k]; }
            }
        }
        cerr << '\n';

        ofstream os(read_args.GetProfilesName().c_str());
        os << "site\tA\tC\tD\tE\tF\tG\tH\tI\tK\tL\tM\tN\tP\tQ\tR\tS\tT\tV\tW\tY\n";
        for (int i = 0; i < model.GetNsite(); i++) {
            os << i + 1;
            for (auto &aa : sitestat[i]) {
                aa /= size;
                os << '\t' << aa;
            }
            os << '\n';
        }
        cerr << "mean site-specific profiles in " << read_args.GetProfilesName() << "\n";
        cerr << '\n';
    } else if (!read_args.chain_omega.getValue().empty()) {
        string chain_omega = read_args.chain_omega.getValue();
        ifstream is_chain_omega{chain_omega + ".param"};
        ChainDriver::fake_read(is_chain_omega);
        AAMutSelMultipleOmegaModel model_omega(is_chain_omega);
        ChainReader cr_omega{model_omega, chain_omega + ".chain"};
        cr_omega.skip(burnin);
        if (not model_omega.FlatFitness()) {
            cerr << "Error: the chain " << chain_omega
                 << " does not have the option `flatfitness`.\n";
            exit(1);
        }
        if (not model.IsDataEqual(model_omega)) {
            cerr << "Error: the two chains are obtained from different dataset.\n";
            exit(1);
        }
        vector<vector<double>> sites_w(model.GetNsite());
        vector<vector<double>> sites_w_knot(model.GetNsite());
        vector<vector<double>> sites_wA_phy(model.GetNsite());
        vector<double> gene_w{};
        vector<double> gene_w_knot{};
        vector<double> gene_wA_phy{};

        for (int step = 0; step < size; step++) {
            cerr << '.';
            cr.skip(every);
            cr_omega.skip(every);
            double tot_w{0.0}, tot_w0{0.0}, tot_wA{0.0};
            for (int site = 0; site < model.GetNsite(); site++) {
                double w = model_omega.GetSiteOmega(site);
                double w0 = model.GetPredictedSiteOmegaKnot(site);
                sites_w[site].push_back(w);
                sites_w_knot[site].push_back(w0);
                sites_wA_phy[site].push_back(w - w0);
                tot_w += w;
                tot_w0 += w0;
                tot_wA += w - w0;
            }
            gene_w.push_back(tot_w / model.GetNsite());
            gene_w_knot.push_back(tot_w0 / model.GetNsite());
            gene_wA_phy.push_back(tot_wA / model.GetNsite());
        }
        cerr << '\n';

        string filename{chain_name + ".omegaA.tsv"};
        ofstream os(filename.c_str());
        os << "site\tω\tω₀\tωᴬ=ω-ω₀\tp(ωᴬ>0)\n";

        double post_w = accumulate(gene_w.begin(), gene_w.end(), 0.0) / size;
        double post_w0 = accumulate(gene_w_knot.begin(), gene_w_knot.end(), 0.0) / size;
        double post_wA = accumulate(gene_wA_phy.begin(), gene_wA_phy.end(), 0.0) / size;
        assert(abs(post_wA - (post_w - post_w0) < 1e-8));
        double proba = static_cast<double>(count_if(gene_wA_phy.begin(), gene_wA_phy.end(),
                           [](double wA) { return wA > 0.0; })) /
                       static_cast<double>(size);
        assert(proba >= 0.0 and proba <= 1.0);
        os << "Mean\t" << post_w << '\t' << post_w0 << '\t' << post_wA << '\t' << proba << '\n';

        for (int i = 0; i < model.GetNsite(); i++) {
            post_w = accumulate(sites_w[i].begin(), sites_w[i].end(), 0.0) / size;
            post_w0 = accumulate(sites_w_knot[i].begin(), sites_w_knot[i].end(), 0.0) / size;
            post_wA = accumulate(sites_wA_phy[i].begin(), sites_wA_phy[i].end(), 0.0) / size;
            assert(abs(post_wA - (post_w - post_w0) < 1e-8));
            proba = static_cast<double>(count_if(sites_wA_phy[i].begin(), sites_wA_phy[i].end(),
                        [](double wA) { return wA > 0.0; })) /
                    static_cast<double>(size);
            assert(proba >= 0.0 and proba <= 1.0);
            os << i + 1 << '\t' << post_w << '\t' << post_w0 << '\t' << post_wA << '\t' << proba
               << '\n';
        }
        cerr << '\n';
    } else if (read_args.omega.getValue() or read_args.omega_knot.getValue()) {
        double ci = stod(read_args.confidence_interval.getValue());
        vector<vector<double>> omega(model.GetNsite());
        vector<double> gene_omega{};
        double upper = max(ci, 1.0 - ci);
        double lower = min(ci, 1.0 - ci);

        for (int step = 0; step < size; step++) {
            cerr << '.';
            cr.skip(every);
            double mean{0.0};
            for (int site = 0; site < model.GetNsite(); site++) {
                double val = read_args.omega_knot.getValue() ? model.GetPredictedSiteOmegaKnot(site)
                                                             : model.GetSiteOmega(site);
                omega[site].push_back(val);
                mean += val;
            }
            gene_omega.push_back(mean / model.GetNsite());
        }
        cerr << '\n';
        string filename{chain_name + ".ci" + read_args.confidence_interval.getValue() + ".tsv"};
        ofstream os(filename.c_str());
        string omega_str = read_args.omega_knot.getValue() ? "ω₀" : "ω";
        string lower_str = "CI[" + to_string(lower) + "]";
        string upper_str = "CI[" + to_string(upper) + "]";
        os << "site\t" + lower_str + "\t" + omega_str + "\t" + upper_str + "\n";

        double mean = accumulate(gene_omega.begin(), gene_omega.end(), 0.0) / size;
        sort(gene_omega.begin(), gene_omega.end());
        auto pt_up = static_cast<size_t>(upper * gene_omega.size());
        if (pt_up >= gene_omega.size()) { pt_up = gene_omega.size() - 1; }
        double up = gene_omega.at(pt_up);
        double down = gene_omega.at(static_cast<size_t>(lower * gene_omega.size()));
        os << "Mean\t" << down << '\t' << mean << '\t' << up << '\n';

        for (int i = 0; i < model.GetNsite(); i++) {
            mean = accumulate(omega[i].begin(), omega[i].end(), 0.0) / size;
            sort(omega[i].begin(), omega[i].end());
            pt_up = static_cast<size_t>(upper * omega[i].size());
            if (pt_up >= omega[i].size()) { pt_up = omega[i].size() - 1; }
            up = omega[i].at(pt_up);
            down = omega[i].at(static_cast<size_t>(lower * omega[i].size()));
            os << i + 1 << '\t' << down << '\t' << mean << '\t' << up << '\n';
        }
        cerr << '\n';
    } else if (read_args.nuc.getValue()) {
        vector<vector<double>> rates(Nnuc * Nnuc);
        for (int step = 0; step < size; step++) {
            cerr << '.';
            cr.skip(every);
            for (int i = 0; i < Nnuc; i++) {
                for (int j = 0; j < Nnuc; j++) {
                    if (i != j) {
                        int r = i * Nnuc + j;
                        rates[r].push_back(model.GetNucRate(i, j));
                    }
                }
            }
        }
        cerr << '\n';
        string filename{chain_name + ".nucmatrix.tsv"};
        ofstream os(filename.c_str());
        os << "Name\tRate\n";
        for (int i = 0; i < Nnuc; i++) {
            for (int j = 0; j < Nnuc; j++) {
                if (i != j) {
                    int r = i * Nnuc + j;
                    double q_mean = accumulate(rates.at(r).begin(), rates.at(r).end(), 0.0) / size;
                    os << "q_" << DNAletters[i] << "_" << DNAletters[j] << "\t" << q_mean << '\n';
                }
            }
        }
        cerr << '\n';
    } else {
        double omega_pp = stod(read_args.omega_pp.getValue());
        vector<double> omegappgto(model.GetNsite(), 0);
        vector<double> omega(model.GetNsite(), 0);

        for (int step = 0; step < size; step++) {
            cerr << '.';
            cr.skip(every);
            for (int site = 0; site < model.GetNsite(); site++) {
                omega[site] += model.GetSiteOmega(site);
                if (model.GetSiteOmega(site) > omega_pp) { omegappgto[site]++; }
            }
        }
        cerr << '\n';

        string filename{chain_name + ".omegappgt" + read_args.omega_pp.getValue() + ".tsv"};
        ofstream os(filename.c_str());
        if (model.FlatFitness()) {
            os << "site\tω\tp(ω>" << read_args.omega_pp.getValue() << ")\n";
        } else {
            os << "site\tω*\tp(ω⁎>" << read_args.omega_pp.getValue() << ")\n";
        }

        for (int i = 0; i < model.GetNsite(); i++) {
            os << i + 1 << '\t' << omega[i] / size << '\t' << omegappgto[i] / size << '\n';
        }
        cerr << "Posterior prob of omega greater than " << read_args.omega_pp.getValue() << " in "
             << filename << "\n";
        cerr << '\n';
    }
}