#include "submodels/globom_model.hpp"

int main(int argc, char* argv[]) {
    // parsing command-line arguments
    ChainCmdLine cmd{argc, argv, "SingleOmega", ' ', "0.1"};
    InferenceAppArgParse args(cmd);
    cmd.parse();

    // input data
    auto data = prepare_data(args.alignment.getValue(), args.treefile.getValue());

    // random generator
    auto gen = make_generator();

    // model
    auto model = globom::make(data, gen);

    // move success stats
    MoveStatsRegistry ms;

    // move schedule
    auto touch_matrices = [&model]() {
        auto& nuc_matrix_proxy = get<nuc_rates, matrix_proxy>(model);
        nuc_matrix_proxy.gather();
        codon_submatrix_(model).SetOmega(get<global_omega, omega, value>(model));
        codon_submatrix_(model).CorruptMatrix();
    };

    auto scheduler = make_move_scheduler([&gen, &touch_matrices, &model, &ms]() {
        // move phyloprocess
        touch_matrices();
        phyloprocess_(model).Move(1.0);

        // move omega
        for (int rep = 0; rep < 30; rep++) {
            // move branch lengths
            bl_suffstats_(model).gather();
            branchlengths_sm::gibbs_resample(branch_lengths_(model), bl_suffstats_(model), gen);

            // move omega
            path_suffstats_(model).gather();
            omegapath_suffstats_(model).gather();
            omega_sm::gibbs_resample(global_omega_(model), omegapath_suffstats_(model), gen);

            // move nuc rates
            touch_matrices();
            nucpath_suffstats_(model).Clear();
            nucpath_suffstats_(model).AddSuffStat(
                codon_submatrix_(model), path_suffstats_(model).get());
            auto nucrates_logprob = [&model]() {
                return nucpath_suffstats_(model).GetLogProb(
                    get<nuc_rates, matrix_proxy>(model).get(), codon_statespace_(model));
            };
            auto touch_nucmatrix = [&model]() { get<nuc_rates, matrix_proxy>(model).gather(); };
            nucrates_sm::move_exch_rates(nuc_rates_(model), {0.1, 0.03, 0.01}, nucrates_logprob,
                touch_nucmatrix, gen, ms("exch_rates"));
            nucrates_sm::move_eq_freqs(nuc_rates_(model), {0.1, 0.03}, nucrates_logprob,
                touch_nucmatrix, gen, ms("eq_freqs"));
            touch_matrices();
        }
    });

    // initializing components
    ChainDriver chain_driver{cmd.chain_name(), args.every.getValue(), args.until.getValue()};

    ConsoleLogger console_logger;
    // ChainCheckpoint chain_checkpoint(cmd.chain_name() + ".param", chain_driver, model);
    StandardTracer trace(model, cmd.chain_name());

    // registering components to chain driver
    chain_driver.add(scheduler);
    chain_driver.add(console_logger);
    // chain_driver.add(chain_checkpoint);
    chain_driver.add(trace);
    chain_driver.add(ms);

    // launching chain!
    chain_driver.go();
}
