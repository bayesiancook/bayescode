/*Copyright or © or Copr. CNRS (2019). Contributors:
- Vincent Lanore. vincent.lanore@gmail.com

This software is a computer program whose purpose is to provide a set of C++ data structures and
functions to perform Bayesian inference with MCMC algorithms.

This software is governed by the CeCILL-C license under French law and abiding by the rules of
distribution of free software. You can use, modify and/ or redistribute the software under the terms
of the CeCILL-C license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info".

As a counterpart to the access to the source code and rights to copy, modify and redistribute
granted by the license, users are provided only with a limited warranty and the software's author,
the holder of the economic rights, and the successive licensors have only limited liability.

In this respect, the user's attention is drawn to the risks associated with loading, using,
modifying and/or developing or reproducing the software by the user in light of its specific status
of free software, that may mean that it is complicated to manipulate, and that also therefore means
that it is reserved for developers and experienced professionals having in-depth computer knowledge.
Users are therefore encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or data to be ensured and,
more generally, to use and operate it in the same conditions as regards security.

The fact that you are presently reading this means that you have had knowledge of the CeCILL-C
license and that you accept its terms.*/

#include <iostream>
#include "basic_moves.hpp"
#include "distributions/exponential.hpp"
#include "distributions/gamma.hpp"
#include "distributions/poisson.hpp"
#include "global/logging.hpp"
#include "mcmc_utils.hpp"
#include "operations/backup.hpp"
#include "operations/logprob.hpp"
#include "operations/raw_value.hpp"
#include "operations/set_value.hpp"
#include "structure/View.hpp"
#include "structure/array_utils.hpp"
#include "suffstat_utils.hpp"
#include "tagged_tuple/src/fancy_syntax.hpp"
using namespace std;

TOKEN(alpha)
TOKEN(mu)
TOKEN(lambda)
TOKEN(K)

auto poisson_gamma(size_t size, size_t size2) {
    auto alpha = make_node<exponential>(1.0);
    auto mu = make_node<exponential>(1.0);
    auto lambda = make_node_array<gamma_ss>(size, n_to_one(alpha), n_to_one(mu));
    auto K = make_node_matrix<poisson>(size, size2,
                                       [& v = get<value>(lambda)](int i, int) { return v[i]; });
    // clang-format off
    return make_model(
         alpha_ = move(alpha),
            mu_ = move(mu),
        lambda_ = move(lambda),
             K_ = move(K)
    );  // clang-format on
}

template <class Node, class MB, class Gen, class... IndexArgs>
void scaling_move(Node& node, MB blanket, Gen& gen, IndexArgs... args) {
    auto index = make_index(args...);
    auto bkp = backup(node, index);
    double logprob_before = logprob(blanket);
    double log_hastings = scale(raw_value(node, index), gen);
    bool accept = decide(logprob(blanket) - logprob_before + log_hastings, gen);
    if (!accept) { restore(node, bkp, index); }
}

int main() {
    auto gen = make_generator();

    constexpr size_t nb_it{100'000}, len_lambda{5}, len_K{6};
    auto m = poisson_gamma(len_lambda, len_K);

    auto v = make_view<alpha, mu, lambda, K>(m);
    draw(v, gen);
    // display node value in stdout
    INFO("Alpha = {}", get<alpha, value>(m));
    INFO("Mu = {}", get<mu, value>(m));
    INFO("Lambda = {}", vector_to_string(get<lambda, value>(m)));
    // INFO("K = {}", get<K, value>(m));

//    set_value(K_(m), {{1, 2, 1}, {1, 2, 2}, {1, 2, 1}, {2, 1, 2}, {2, 1, 3}});

    double alpha_sum{0}, mu_sum{0};
    vector<double> lambda_sum (len_lambda, 0.0);

    for (size_t it = 0; it < nb_it; it++) {
        scaling_move(alpha_(m), make_view<alpha, lambda>(m), gen);
        scaling_move(mu_(m), make_view<mu, lambda>(m), gen);
        alpha_sum += raw_value(alpha_(m));
        mu_sum += raw_value(mu_(m));

        for (size_t i = 0; i < len_lambda; i++) {
            auto lambda_mb = make_view(make_ref<K>(m, i), make_ref<lambda>(m, i));
            scaling_move(lambda_(m), lambda_mb, gen, i);
            lambda_sum[i] += raw_value(lambda_(m), i);
        }
    }
    INFO("RESULTS OF THE MCMC: ");
    INFO("Alpha = {}", alpha_sum / float(nb_it));
    INFO("Mu = {}", mu_sum / float(nb_it));

    for (size_t i = 0; i < len_lambda; i++) {
      lambda_sum[i] = lambda_sum[i]/float(nb_it) ;
    }

    //std::cout << "Mean lambda = " << lambda_sum / (float(nb_it) * len_lambda) << std::endl;
    INFO("Lambda = {}", vector_to_string(lambda_sum));

}
