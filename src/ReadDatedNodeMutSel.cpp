#include <cmath>
#include <fstream>
#include "DatedNodeMutSelModel.hpp"
#include "components/ChainDriver.hpp"
#include "components/ChainReader.hpp"
#include "components/ReadArgParse.hpp"
#include "components/stats_posterior.hpp"
#include "tclap/CmdLine.h"
#include "tree/export.hpp"

using namespace std;
using namespace TCLAP;


class ReadNodeMutSelArgParse : public ReadArgParse {
  public:
    explicit ReadNodeMutSelArgParse(CmdLine &cmd) : ReadArgParse(cmd) {}

    SwitchArg ss{"s", "ss",
        "Computes the mean posterior site-specific amino-acid equilibrium frequencies"
        "(amino-acid fitness profiles). "
        "Results are written in {chain_name}.siteprofiles by default (optionally use the --output argument "
        " to specify a different output path).",
        cmd};
    SwitchArg newick{"t", "newick",
        "Computes the mean posterior node-specific entries of the multivariate Brownian process. "
        "Each entry of the multivariate Brownian process is written in a newick extended (.nhx) "
        "format file."
        "For each trait, results are written in {chain_name}.{trait}.nhx by default (optionally "
        "use the --output argument to specify a different output path).",
        cmd};
    SwitchArg cov{"c", "cov",
        "Computes the mean posterior covariance matrix, precision matrix and correlation matrix. "
        "Results are written in {chain_name}.cov by default (optionally use the --output argument "
        "to specify a different output path).",
        cmd};
};


int main(int argc, char *argv[]) {
    CmdLine cmd{"DatedMutSel", ' ', "0.1"};
    ReadNodeMutSelArgParse read_args(cmd);
    cmd.parse(argc, argv);

    string chain_name = read_args.GetChainName();
    int burnin = read_args.GetBurnIn();
    int every = read_args.GetEvery();
    int size = read_args.GetSize();

    ifstream is{chain_name + ".param"};
    unique_ptr<DatedNodeMutSelModel> model = nullptr;
    new ChainDriver(is);
    is >> model;
    ChainReader cr(*model, chain_name + ".chain");

    cr.skip(burnin);
    cerr << size << " points to read\n";

    if (read_args.GetPpred()) {
        for (int i = 0; i < size; i++) {
            cerr << '.';
            cr.skip(every);
            model->PostPred("ppred_" + chain_name + "_" + to_string(i) + ".ali");
        }
        cerr << '\n';
    } if (read_args.trace.getValue()) {
        string file_name = read_args.OutputFile(".trace.tsv");
        recompute_trace<DatedNodeMutSelModel>(*model, cr, file_name, every, size);
    } else if (read_args.ss.getValue()) {
        string file_name = read_args.OutputFile(".siteprofiles");
        vector<vector<double>> sitestat(model->GetNsite(), {0});

        for (int step = 0; step < size; step++) {
            cerr << '.';
            cr.skip(every);
            for (int i = 0; i < model->GetNsite(); i++) {
                vector<double> const &profile = model->GetProfile(i);
                if (sitestat[i].size() != profile.size()) { sitestat[i].resize(profile.size(), 0); }
                for (unsigned k{0}; k < profile.size(); k++) { sitestat[i][k] += profile[k]; }
            }
        }
        cerr << '\n';

        ofstream os(file_name);
        os << "site\tA\tC\tD\tE\tF\tG\tH\tI\tK\tL\tM\tN\tP\tQ\tR\tS\tT\tV\tW\tY\n";
        for (int i = 0; i < model->GetNsite(); i++) {
            os << i + 1;
            for (auto &aa : sitestat[i]) {
                aa /= size;
                os << '\t' << aa;
            }
            os << '\n';
        }
        cerr << "mean site-specific profiles in " << file_name << "\n";
        cerr << '\n';
    } else if (read_args.cov.getValue()) {
        string file_name = read_args.OutputFile(".cov");
        ofstream os(file_name);
        os << "entries are in the following order:" << endl;

        for (int dim = 0; dim < model->GetDimension(); dim++) {
            os << model->GetDimensionName(dim) << endl;
        }

        EMatrix posterior_prob = EMatrix::Zero(model->GetDimension(), model->GetDimension());
        EMatrix cov_matrix = EMatrix::Zero(model->GetDimension(), model->GetDimension());
        for (int step = 0; step < size; step++) {
            cerr << '.';
            cr.skip(every);
            EMatrix cov_matrix_chain = model->GetPrecisionMatrix().inverse();
            cov_matrix += cov_matrix_chain;
            for (int i = 0; i < model->GetDimension(); i++) {
                for (int j = 0; j < model->GetDimension(); j++) {
                    if (cov_matrix_chain.coeffRef(i, j) > 0) { posterior_prob.coeffRef(i, j) += 1; }
                }
            }
        }
        posterior_prob /= size;
        cov_matrix /= size;

        export_matrix(os, model->GetDimension(), cov_matrix, "covariances");
        EMatrix cor_matrix = cov_matrix;
        for (int i = 0; i < model->GetDimension(); i++) {
            for (int j = 0; j < model->GetDimension(); j++) {
                cor_matrix.coeffRef(i, j) =
                    cov_matrix.coeffRef(i, j) /
                    sqrt(cov_matrix.coeffRef(i, i) * cov_matrix.coeffRef(j, j));
            }
        }
        export_matrix(os, model->GetDimension(), cor_matrix, "correlation coefficients");
        export_matrix(os, model->GetDimension(), posterior_prob, "posterior probs", false);

        cerr << endl << "matrices in " << file_name << "." << endl;
    } else if (read_args.newick.getValue()) {
        vector<vector<vector<double>>> dim_node_traces(model->GetDimension());
        vector<vector<double>> branch_times(model->GetTree().nb_nodes());
        vector<vector<double>> branch_length(model->GetTree().nb_nodes());
        vector<vector<double>> leaves_theta(model->GetTree().nb_nodes());
        vector<vector<double>> contrast_pop_size(model->GetTree().nb_nodes());

        for (int dim{0}; dim < model->GetDimension(); dim++) {
            dim_node_traces[dim].resize(model->GetTree().nb_nodes());
        }

        for (int step = 0; step < size; step++) {
            cerr << '.';
            cr.skip(every);

            model->UpdateBranches(true);
            for (Tree::NodeIndex node = 0; node < Tree::NodeIndex(model->GetTree().nb_nodes());
                 node++) {
                if (!model->GetTree().is_root(node)) {
                    double branch_time = model->GetBranchTime(node);
                    assert(branch_time >= 0);
                    branch_times[node].push_back(branch_time);
                    branch_length[node].push_back(model->GetBranchLength(node));
                    contrast_pop_size[node].push_back(model->GetContrast(node, dim_pop_size));
                }
                if (model->PolymorphismAware() and model->GetTree().is_leaf(node)) {
                    leaves_theta[node].push_back(model->GetTheta(node));
                }
                for (int dim{0}; dim < model->GetDimension(); dim++) {
                    dim_node_traces[dim][node].push_back(model->GetExpBrownianEntry(node, dim));
                }
            }
        }
        cerr << '\n';

        ExportTree base_export_tree(model->GetTree());
        for (Tree::NodeIndex node = 0; node < Tree::NodeIndex(model->GetTree().nb_nodes());
             node++) {
            if (!model->GetTree().is_root(node)) {
                base_export_tree.set_tag(node, "length", to_string(mean(branch_times[node])));
            }
        }

        export_tree(
            base_export_tree, "ContrastPopulationSize", read_args.OutputFile(), contrast_pop_size);
        export_tree(base_export_tree, "BranchLength", read_args.OutputFile(), branch_length);
        export_tree(base_export_tree, "BranchTime", read_args.OutputFile(), branch_times);
        if (model->PolymorphismAware()) {
            export_tree(base_export_tree, "Theta", read_args.OutputFile(), leaves_theta);
        }
        for (int dim{0}; dim < model->GetDimension(); dim++) {
            export_tree(base_export_tree, model->GetDimensionName(dim), read_args.OutputFile(),
                dim_node_traces[dim]);
        }
    } else {
        stats_posterior<DatedNodeMutSelModel>(*model, cr, every, size);
    }
}