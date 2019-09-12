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

template <class Gen>
auto make_nuc_rates(const std::vector<double>& nucrelratehypercenter, double nucrelratehyperinvconc,
    const std::vector<double>& nucstathypercenter, double nucstathyperinvconc, Gen& gen) {
    /* -- */
    auto exchangeability_rates = make_vector_node<dirichlet_cic>(
        6, std::move(nucrelratehypercenter), std::move(nucrelratehyperinvconc));
    draw(exchangeability_rates, gen);

    auto equilibrium_frequencies = make_vector_node<dirichlet_cic>(
        4, std::move(nucstathypercenter), std::move(nucstathyperinvconc));
    draw(equilibrium_frequencies, gen);

    auto nuc_matrix = std::make_unique<GTRSubMatrix>(
        4, get<value>(exchangeability_rates), get<value>(equilibrium_frequencies), true);

    return make_model(                                   //
        exch_rates_ = std::move(exchangeability_rates),  //
        eq_freq_ = std::move(equilibrium_frequencies),   //
        nuc_matrix_ = std::move(nuc_matrix));
}