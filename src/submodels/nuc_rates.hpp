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
        6, std::move(nucrelratecenter), std::move(nucrelrateinvconc));
    draw(exchangeability_rates, gen);
    DEBUG("GTR model: exchangeability rates are {}.",
        vector_to_string(get<value>(exchangeability_rates)));

    auto equilibrium_frequencies =
        make_vector_node<dirichlet_cic>(4, std::move(nucstatcenter), std::move(nucstatinvconc));
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