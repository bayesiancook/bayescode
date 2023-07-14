#include <cmath>
#include <fstream>
#include "DatedNodeModel.hpp"
#include "components/ChainDriver.hpp"
#include "components/ChainReader.hpp"
#include "components/ReadArgParse.hpp"
#include "components/csv_parser.hpp"
#include "components/stats_posterior.hpp"
#include "tclap/CmdLine.h"
#include "tree/export.hpp"

using namespace std;
using namespace TCLAP;


class ReadNodeArgParse : public ReadArgParse {
  public:
    explicit ReadNodeArgParse(CmdLine& cmd) : ReadArgParse(cmd) {}

    SwitchArg newick{"t", "newick",
        "Computes the mean posterior node-specific entries of the multivariate Brownian process. "
        "Each entry of the multivariate Brownian process is written in a newick extended (.nhx) "
        "format file.",
        cmd};

    SwitchArg cov{"c", "cov", "Computes the mean posterior covariance matrix.", cmd};

    ValueArg<string> var_within{"", "var_within",
        "A .tsv file containing the within-population variance of each trait.", false, "", "string",
        cmd};

    ValueArg<string> output{"o", "output", "Output file name (optional)", false, "", "string", cmd};

    string OutputFile(const string& default_suffix = "") {
        if (!output.getValue().empty()) {
            return output.getValue();
        } else {
            cout << "No output file name specified, using default: "
                 << GetChainName() + default_suffix << endl;
            return GetChainName() + default_suffix;
        }
    }
};

class VarWithinSample {
  public:
    VarWithinSample() = default;

    VarWithinSample(const CSVParser& var_within, const string& trait_name) {
        size_t taxa_col_index = var_within.find_column_index("TaxonName");
        if (taxa_col_index == string::npos) {
            cerr << "Error: Could not find column 'TaxonName'" << endl;
            exit(1);
        }
        size_t ps_index = var_within.find_column_index("Nucleotide_diversity");
        if (ps_index == string::npos) {
            cerr << "Error: Could not find column 'Nucleotide_diversity'" << endl;
            exit(1);
        }
        size_t trait_var_index = var_within.find_column_index(trait_name + "_variance");
        if (trait_var_index == string::npos) {
            cerr << "Error: Could not find column '" << trait_name << "_variance'" << endl;
            exit(1);
        }

        size_t h_index = var_within.find_column_index(trait_name + "_heritability");
        size_t h_low_index = var_within.find_column_index(trait_name + "_heritability_lower");
        size_t h_up_index = var_within.find_column_index(trait_name + "_heritability_upper");

        random_sampling = (h_low_index != string::npos && h_up_index != string::npos);
        if (random_sampling) {
            cout << "Random sampling of within-population variance using heritability bounds"
                    " for trait '"
                 << trait_name << "'" << endl;
        } else if (h_index != string::npos) {
            cout << "Using heritability for trait '" << trait_name << "'" << endl;
        } else {
            cout << "No heritability information for trait '" << trait_name << "'" << endl;
        }
        double sum_var_within = 0.0;
        for (size_t taxa = 0; taxa < var_within.nbr_rows(); ++taxa) {
            double taxa_var = var_within.get_value_as_double(trait_var_index, taxa);
            if (isnan(taxa_var)) {
                continue;
            } else if (taxa_var < 0.0) {
                cerr << "Error: Within-population variance is negative (" << taxa_var
                     << ") for taxa " << var_within.get_value(taxa_col_index, taxa) << "." << endl;
                exit(1);
            }
            double taxa_nuc_div = var_within.get_value_as_double(ps_index, taxa);
            if (isnan(taxa_nuc_div)) {
                continue;
            } else if (taxa_nuc_div < 0.0) {
                cerr << "Error: 'Nucleotide_diversity' is negative (" << taxa_nuc_div
                     << ") for taxa " << var_within.get_value(taxa_col_index, taxa) << "." << endl;
                exit(1);
            }
            double r = taxa_var / taxa_nuc_div;
            if (random_sampling) {
                double taxa_h_low = var_within.get_value_as_double(h_low_index, taxa);
                double taxa_h_up = var_within.get_value_as_double(h_up_index, taxa);
                if (taxa_h_low < 0.0 or taxa_h_up > 1.0) {
                    cerr << "Error: Heritability bounds are not in [0, 1] for taxa "
                         << var_within.get_value(taxa_col_index, taxa) << ", skipping it." << endl;
                    continue;
                }
                var_within_lower_bound.push_back(r * taxa_h_low);
                var_within_upper_bound.push_back(r * taxa_h_up);
                sum_var_within += r * (taxa_h_low + taxa_h_up) / 2.0;
            } else if (h_index != string::npos) {
                double h = var_within.get_value_as_double(h_index, taxa);
                if (isnan(h) or h < 0.0 or h > 1.0) {
                    cerr << "Error: Heritability is not in [0, 1] for taxa "
                         << var_within.get_value(taxa_col_index, taxa) << ", skipping taxa."
                         << endl;
                    continue;
                }
                sum_var_within += r * h;
            } else {
                sum_var_within += r;
            }
            n_taxa++;
        }
        mean_var_within = sum_var_within / n_taxa;
        cout << "Mean within-population variance: " << mean_var_within << " for trait "
             << trait_name << endl;
    }

    double GetSample() {
        if (random_sampling) {
            double sum_var = 0.0;
            for (size_t i = 0; i < var_within_lower_bound.size(); ++i) {
                double var_within_low = var_within_lower_bound[i];
                double var_within_high = var_within_upper_bound[i];
                double rand = Random::Uniform();
                if (rand < 0.0 || rand > 1.0) {
                    cerr << "Error: Random number is not between 0 and 1 (" << rand << ")" << endl;
                    exit(1);
                }
                sum_var += rand * (var_within_high - var_within_low) + var_within_low;
            }
            return sum_var / n_taxa;
        } else {
            return mean_var_within;
        }
    }

    int GetNbrTaxa() const { return n_taxa; }

  private:
    bool random_sampling{false};
    double mean_var_within{0.0};
    int n_taxa{0};
    vector<double> var_within_lower_bound;
    vector<double> var_within_upper_bound;
};

tuple<double, double> quantile(vector<double>& sorted_v, double q) {
    if (sorted_v.empty()) { return make_tuple(0.0, 0.0); }
    if (q < 0.0 or q > 1.0) {
        cerr << "quantile is not in [0,1]\n";
        exit(1);
    }
    if (q > 0.5) { q = 1.0 - q; }
    auto lo = static_cast<size_t>(floor(q * static_cast<double>(sorted_v.size())));
    auto hi = static_cast<size_t>(ceil((1.0 - q) * static_cast<double>(sorted_v.size()))) - 1;
    assert(lo <= hi);
    assert(hi < sorted_v.size());
    return make_tuple(sorted_v.at(lo), sorted_v.at(hi));
}

string replace_last_occurrence(string str, const string& from, const string& to) {
    size_t start_pos = str.rfind(from);
    if (start_pos == string::npos) { return str; }
    str.replace(start_pos, from.length(), to);
    return str;
}

string format_dim(string dim) { return replace_last_occurrence(std::move(dim), "_mean", ""); }

int main(int argc, char* argv[]) {
    CmdLine cmd{"DatedMutSel", ' ', "0.1"};
    ReadNodeArgParse read_args(cmd);
    cmd.parse(argc, argv);

    string chain_name = read_args.GetChainName();
    int burnin = read_args.GetBurnIn();
    int every = read_args.GetEvery();
    int size = read_args.GetSize();

    ifstream is{chain_name + ".param"};
    ChainDriver* fake_read = nullptr;
    unique_ptr<DatedNodeModel> model = nullptr;
    fake_read = new ChainDriver(is);
    is >> model;
    ChainReader cr(*model, chain_name + ".chain");

    cr.skip(burnin);
    cerr << size << " points to read\n";

    if (read_args.trace.getValue()) {
        recompute_trace<DatedNodeModel>(*model, cr, chain_name, every, size);
    } else if (not read_args.var_within.getValue().empty()) {
        CSVParser tsv_war_within(read_args.var_within.getValue(), '\t');
        unordered_map<string, VarWithinSample> var_within_dict;

        for (int dim = 0; dim < model->GetDimension(); dim++) {
            string dim_name = format_dim(model->GetDimensionName(dim));
            var_within_dict[dim_name] = VarWithinSample(tsv_war_within, dim_name);
        }

        EVector var_between_trace = EVector::Zero(model->GetDimension());
        EVector var_within_trace = EVector::Zero(model->GetDimension());
        EVector ratio_trace = EVector::Zero(model->GetDimension());
        EVector posterior_prob = EVector::Zero(model->GetDimension());
        vector<vector<double>> ratio_trace_as_vector =
            vector<vector<double>>(model->GetDimension());
        for (int step = 0; step < size; step++) {
            cerr << '.';
            cr.skip(every);
            EVector var_between_chain = model->GetCovarianceMatrix().diagonal() / 4.0;
            for (int dim = 0; dim < model->GetDimension(); dim++) {
                string name = format_dim(model->GetDimensionName(dim));
                if (var_between_chain(dim) < 0) {
                    cerr << "error: negative variance\n";
                    exit(1);
                }

                double var_within_sample = var_within_dict.at(name).GetSample();
                var_within_trace(dim) += var_within_sample;
                var_between_trace(dim) += var_between_chain(dim);
                double ratio = var_between_chain(dim) / var_within_sample;
                ratio_trace(dim) += ratio;
                ratio_trace_as_vector[dim].push_back(ratio);
                if (ratio > 1.0) { posterior_prob(dim) += 1; }
            }
        }
        posterior_prob /= size;
        var_between_trace /= size;
        var_within_trace /= size;
        ratio_trace /= size;

        cerr << endl;
        string file_name = read_args.OutputFile(".ratio.tsv");
        ofstream os(file_name);
        os << "trait\tnbr_taxa_between\tvar_between\tnbr_taxa_within\tvar_within\tratio\t"
              "pp_ratio_greater_1\tpp_ratio_lower_1\tratio_CI95\tratio_CI99\n";
        for (int dim = 0; dim < model->GetDimension(); dim++) {
            string dim_name = format_dim(model->GetDimensionName(dim));
            vector<double> ratio_trace_sorted = ratio_trace_as_vector[dim];
            sort(ratio_trace_sorted.begin(), ratio_trace_sorted.end());
            auto ratio_CI95 = quantile(ratio_trace_sorted, 0.025);
            auto ratio_CI99 = quantile(ratio_trace_sorted, 0.005);

            os << dim_name << '\t' << model->GetNtaxa() << '\t' << var_between_trace(dim) << '\t'
               << var_within_dict.at(dim_name).GetNbrTaxa() << '\t' << var_within_trace(dim) << '\t'
               << ratio_trace(dim) << '\t' << posterior_prob(dim) << '\t'
               << 1 - posterior_prob(dim);
            os << '\t' << get<0>(ratio_CI95) << "--" << get<1>(ratio_CI95);
            os << '\t' << get<0>(ratio_CI99) << "--" << get<1>(ratio_CI99) << '\n';

            cout << "Trait " << format_dim(model->GetDimensionName(dim)) << ":" << endl;
            cout << " - var_between = " << var_between_trace(dim) << endl;
            cout << " - var_within = " << var_within_trace(dim) << endl;
            cout << " - ρ = var_between / var_within = " << ratio_trace(dim) << endl;
            cout << " - posterior P[ρ > 1]  = " << posterior_prob(dim) << endl;
            cout << " - posterior P[ρ <= 1] = " << 1 - posterior_prob(dim) << endl;
        }
        os.close();
        cout << "Results written to " << file_name << "." << endl;
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
            EMatrix cov_matrix_chain = model->GetCovarianceMatrix();
            for (int i = 0; i < model->GetDimension(); i++) {
                if (cov_matrix_chain(i, i) < 0) {
                    std::cerr << "error: negative variance\n";
                    exit(1);
                }
            }
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

        for (int dim{0}; dim < model->GetDimension(); dim++) {
            dim_node_traces[dim].resize(model->GetTree().nb_nodes());
        }

        for (int step = 0; step < size; step++) {
            cerr << '.';
            cr.skip(every);

            model->Update();
            for (Tree::NodeIndex node = 0; node < Tree::NodeIndex(model->GetTree().nb_nodes());
                 node++) {
                if (!model->GetTree().is_root(node)) {
                    double branch_time = model->GetBranchTime(node);
                    assert(branch_time >= 0);
                    branch_times[node].push_back(branch_time);
                }
                for (int dim{0}; dim < model->GetDimension(); dim++) {
                    dim_node_traces[dim][node].push_back(model->GetBrownianEntry(node, dim));
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
        for (int dim{0}; dim < model->GetDimension(); dim++) {
            export_tree(base_export_tree, model->GetDimensionName(dim), read_args.OutputFile(),
                dim_node_traces[dim]);
        }
    } else {
        stats_posterior<DatedNodeModel>(*model, cr, every, size);
    }
}