#pragma once

#include "bayes_toolbox/src/basic_moves.hpp"
#include "bayes_toolbox/src/distributions/exponential.hpp"
#include "bayes_toolbox/src/distributions/gamma.hpp"
#include "bayes_toolbox/src/operations/draw.hpp"
#include "bayes_toolbox/src/structure/array_utils.hpp"
#include "bayes_toolbox/src/structure/model.hpp"
#include "bayes_toolbox/src/structure/node.hpp"
#include "bayes_toolbox/utils/tagged_tuple/src/fancy_syntax.hpp"
#include "tree/implem.hpp"

TOKEN(omega)

template <class Target, class LogprobLambda, class Gen>
void scaling_move(Target& target, LogprobLambda& logprob, Gen& gen) {
    auto bkp = backup(target);
    double logprob_before = logprob();
    double log_hastings = scale(raw_value(target), gen);
    bool accept = decide(logprob() - logprob_before + log_hastings, gen);
    if (!accept) { restore(target, bkp); }
}

template <class Gen>
auto make_fixed_globom(double mean, double invshape, Gen& gen) {
    auto omega = make_node<gamma_ss>(1. / invshape, mean * invshape);
    draw(omega, gen);
    return make_model(omega_ = std::move(omega));
}

template <class GlobomModel, class Lambda, class Gen>
void move_globom(GlobomModel& model, Lambda& children_logprob, Gen& gen) {
    auto& omega = omega_(model);
    auto full_logprob = [&omega, &children_logprob]() {
        return logprob(omega) + children_logprob();
    };
    scaling_move(omega, full_logprob);
}

TOKEN(bl_array)

/* Array of branch lengths, gamma iid with fixed mean and invshape. */
auto make_branchlength_array(TreeParser& parser, double mean, double invshape) {
    DEBUG("Getting branch lengths from tree");
    const size_t nb_branches = parser.get_tree().nb_nodes();
    auto initial_bl = branch_container_from_parser<double>(
        parser, [](int i, const auto& tree) { return stod(tree.tag(i, "length")); });

    DEBUG("Creating branch length array of gamma nodes (length {})", nb_branches);
    auto bl_array = make_node_array<gamma_ss>(
        nb_branches, n_to_constant(1. / invshape), n_to_constant(mean * invshape));
    set_value(bl_array, initial_bl);

    // return model
    return make_model(bl_array_ = std::move(bl_array));
}