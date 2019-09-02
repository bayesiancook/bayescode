#pragma once

#include "bayes_toolbox/src/distributions/exponential.hpp"
#include "bayes_toolbox/src/distributions/gamma.hpp"
#include "bayes_toolbox/src/structure/model.hpp"
#include "bayes_toolbox/src/structure/node.hpp"
#include "bayes_toolbox/utils/tagged_tuple/src/fancy_syntax.hpp"

struct fixed_param {};
struct variable_param {};
struct shared_param : variable_param {};
struct shrunken_param : variable_param {};

TOKEN(omega)

auto make_omega(fixed_param, double mean, double invshape) {
    auto omega = make_node<gamma_sr>(mean, invshape);  //@todo: wrong gamma

    return make_model(omega_ = std::move(omega));
}

auto make_omega(variable_param, double mean, double invshape) {
    auto omega_hypermean = make_node<exponential>(mean);
    auto omega_hyperinvshape = make_node<exponential>(invshape);
    auto omega = make_node<gamma_sr>(omega_hypermean, omega_hyperinvshape);  //@todo: wrong gamma

    return make_model(omega_ = std::move(omega));
}
