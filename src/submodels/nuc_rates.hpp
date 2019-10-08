#pragma once

#include "bayes_toolbox/src/basic_moves.hpp"
#include "bayes_toolbox/src/distributions/dirichlet.hpp"
#include "bayes_toolbox/src/distributions/exponential.hpp"
#include "bayes_toolbox/src/distributions/gamma.hpp"
#include "bayes_toolbox/src/operations/draw.hpp"
#include "bayes_toolbox/src/operations/logprob.hpp"
#include "bayes_toolbox/src/structure/array_utils.hpp"
#include "bayes_toolbox/src/structure/model.hpp"
#include "bayes_toolbox/src/structure/node.hpp"
#include "bayes_toolbox/utils/tagged_tuple/src/fancy_syntax.hpp"
#include "global/logging.hpp"
#include "lib/GTRSubMatrix.hpp"
#include "tree/implem.hpp"

TOKEN(eq_freq)
TOKEN(exch_rates)
TOKEN(nuc_matrix)

// @todo: move elsewhere
std::vector<double> normalize(const std::vector<double>& vec) {
    double sum = 0;
    for (auto e : vec) { sum += e; }
    std::vector<double> result(vec.size());
    for (size_t i = 0; i < vec.size(); i++) { result[i] = vec[i] / sum; }
    return result;
}


template <class Gen>
auto make_nucleotide_rate(const std::vector<double>& nucrelratecenter, double nucrelrateinvconc,
    const std::vector<double>& nucstatcenter, double nucstatinvconc, Gen& gen) {
    /* -- */
    auto exchangeability_rates = make_vector_node<dirichlet_cic>(
        6, [v = nucrelratecenter]() { return v; }, [v = nucrelrateinvconc]() { return v; });
    draw(exchangeability_rates, gen);
    DEBUG("GTR model: exchangeability rates are {}.",
        vector_to_string(get<value>(exchangeability_rates)));

    auto equilibrium_frequencies = make_vector_node<dirichlet_cic>(
        4, [v = nucstatcenter]() { return v; }, [v = nucstatinvconc]() { return v; });
    draw(equilibrium_frequencies, gen);
    DEBUG("GTR model: equilibrium frequencies are {}.",
        vector_to_string(get<value>(equilibrium_frequencies)));

    auto nuc_matrix = std::make_unique<GTRSubMatrix>(
        4, get<value>(exchangeability_rates), get<value>(equilibrium_frequencies), true);

    return make_model(                                   //
        exch_rates_ = std::move(exchangeability_rates),  //
        eq_freq_ = std::move(equilibrium_frequencies),   //
        nuc_matrix_ = std::move(nuc_matrix));
}

template <class SubModel, class LogProb, class Gen>
auto move_exch_rates(SubModel& model, double tuning, LogProb logprob_children, Gen& gen) {
    auto& target = exch_rates_(model);
    static_assert(is_node_array<std::decay_t<decltype(target)>>::value, "");

    auto bkp = backup(target);
    // DEBUG(
    //     "=====================================================\nMove exch rate: "
    //     "logprob_children={}, logprob={}\n\ttarget={}",
    //     logprob_children(), logprob(target), vector_to_string(get<value>(target)));
    double logprob_before = logprob_children() + logprob(target);
    double log_hastings = profile_move(get<value>(target), tuning, gen);
    double logprob_after = logprob_children() + logprob(target);
    bool accept = decide(logprob_after - logprob_before + log_hastings, gen);
    if (!accept) { restore(target, bkp); }

    // DEBUG("Logprob diff = {}, (before:{}, after:{})\n\t(children_after:{})",
    //     logprob_after - logprob_before, logprob_before, logprob_after, logprob_children());

    // DEBUG(
    //     "Move exch rate(after): logprob_children={}, logprob={}\n\ttarget={}, "
    //     "accept={}\n=====================================================\n",
    //     logprob_children(), logprob(target), vector_to_string(get<value>(target)), accept);
    return static_cast<double>(accept);
}

template <class SubModel, class LogProb, class Gen>
auto move_eq_freqs(SubModel& model, double tuning, LogProb logprob_children, Gen& gen) {
    auto& target = eq_freq_(model);
    static_assert(is_node_array<std::decay_t<decltype(target)>>::value, "");

    auto bkp = backup(target);
    double logprob_before = logprob_children() + logprob(target);
    double log_hastings = profile_move(get<value>(target), tuning, gen);
    double logprob_after = logprob_children() + logprob(target);
    bool accept = decide(logprob_after - logprob_before + log_hastings, gen);
    if (!accept) { restore(target, bkp); }
    return static_cast<double>(accept);
}