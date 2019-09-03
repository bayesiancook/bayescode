#pragma once

#include "bayes_toolbox/src/basic_moves.hpp"
#include "bayes_toolbox/src/distributions/exponential.hpp"
#include "bayes_toolbox/src/distributions/gamma.hpp"
#include "bayes_toolbox/src/operations/draw.hpp"
#include "bayes_toolbox/src/structure/array_utils.hpp"
#include "bayes_toolbox/src/structure/model.hpp"
#include "bayes_toolbox/src/structure/node.hpp"
#include "bayes_toolbox/utils/tagged_tuple/src/fancy_syntax.hpp"
#include "logging.hpp"
#include "tree/implem.hpp"

TOKEN(omega)

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