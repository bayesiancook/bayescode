#pragma once

#include "bayes_toolbox/src/basic_moves.hpp"
#include "bayes_toolbox/src/distributions/exponential.hpp"
#include "bayes_toolbox/src/distributions/gamma.hpp"
#include "bayes_toolbox/src/operations/draw.hpp"
#include "bayes_toolbox/src/operations/logprob.hpp"
#include "bayes_toolbox/src/structure/array_utils.hpp"
#include "bayes_toolbox/src/structure/model.hpp"
#include "bayes_toolbox/src/structure/node.hpp"
#include "bayes_toolbox/utils/tagged_tuple/src/fancy_syntax.hpp"
#include "global/logging.hpp"
#include "tree/implem.hpp"

TOKEN(omega)

struct globom {
    template <class Gen>
    static auto make_fixed(double mean, double invshape, Gen& gen) {
        DEBUG("Making global omega with fixed parameters mean={} and invshape={}", mean, invshape);
        auto omega = make_node<gamma_ss>(1. / invshape, mean * invshape);
        draw(omega, gen);
        return make_model(omega_ = std::move(omega));
    }

    template <class GlobomModel, class Lambda, class Gen>
    static void move(GlobomModel& model, Lambda children_logprob, Gen& gen) {
        auto& omega = omega_(model);
        auto full_logprob = [&omega, &children_logprob]() {
            return logprob(omega) + children_logprob();
        };
        scaling_move(omega, full_logprob, gen);
        DEBUG("Omega = {}", raw_value(omega));
    }
};