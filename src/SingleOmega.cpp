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
    auto scheduler = make_move_scheduler([&gen, &model, &ms]() {
        // move phyloprocess
        globom::touch_matrices(model);
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
            globom::move_nucrates(model, gen, ms);
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
