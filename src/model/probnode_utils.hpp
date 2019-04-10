#pragma once

#include "global/Random.hpp"
#include "global/defs.hpp"

/*--------------------------------------------------------------------------------------------------
  Type aliases to make declarations more informative */

using indicator_t = int;  // FIXME? any benefit to change to char?

template <class T>
using per_cond = vector<T>;

template <class T>
using per_site = vector<T>;

template <class T>
using per_aa = vector<T>;

/*--------------------------------------------------------------------------------------------------
  Hyperparameter structs */

struct hyper_mean_invconc {
    double mean;
    double invconc;
};

istream& operator>>(istream& is, hyper_mean_invconc& hyper) {
    return is >> hyper.mean >> hyper.invconc;
}

struct hyper_mean_invshape {  // FIXME? different from hyper_mean_invconc?
    double mean;
    double invshape;
};

istream& operator>>(istream& is, hyper_mean_invshape& hyper) {
    return is >> hyper.mean >> hyper.invshape;
}

struct hyper_rate {
    double rate;
};

istream& operator>>(istream& is, hyper_rate& hyper) { return is >> hyper.rate; }

struct no_hyper {};

istream& operator>>(istream& is, no_hyper& hyper) { return is; }

/*--------------------------------------------------------------------------------------------------
  Functions to set vectors to specific values */
template <class T>
void set_all_to(vector<T>& data, T value) {
    for (auto&& elem : data) { elem = value; }
}

template <class T>
void set_all_to(vector<vector<T>>& data, T value) {
    for (auto&& subvector : data) { set_all_to(subvector, value); }
}

template <class T>
void set_all_to(vector<T>& data, size_t count, T value) {
    assert(count >= 0);
    if (data.size() > 0) { WARNING("Erasing non-empty vector"); }
    data = vector<T>(count, value);
}

template <class T>
void set_all_to(vector<vector<T>>& data, size_t count_x, size_t count_y, T value) {
    assert(count_x >= 0 and count_y >= 0);
    if (data.size() > 0) { WARNING("Erasing non-empty vector"); }
    data = vector<vector<T>>(count_x, vector<T>(count_y, value));
}

/*--------------------------------------------------------------------------------------------------
  Counting bools set to true in vectors */
size_t count_indicators(vector<indicator_t>& data) {
    int result = 0;
    for (auto ind : data) {
        assert(ind == 0 or ind == 1);
        result += ind > 0 ? 1 : 0;
    }
    return result;
}

size_t count_indicators(vector<vector<indicator_t>>& data) {
    int result = 0;
    for (auto&& subvector : data) { result += count_indicators(subvector); }
    return result;
}

/*--------------------------------------------------------------------------------------------------
  Drawing vector of indicators from iid bernoullis */
void draw_bernoulli_iid(vector<indicator_t>& data, double prob) {
    assert(prob >= 0 and prob <= 1);
    std::bernoulli_distribution distrib(prob);
    for (auto&& elem : data) { elem = distrib(Random::global_gen); }
}

void draw_bernoulli_iid(vector<vector<indicator_t>>& data, const vector<double>& probs) {
    assert(probs.size() == data.size());
    for (size_t i = 0; i < data.size(); i++) { draw_bernoulli_iid(data.at(i), probs.at(i)); }
}

void draw_bernoulli_iid(vector<vector<indicator_t>>& data, double prob) {
    assert(prob >= 0 and prob <= 1);
    for (auto&& subvector : data) { draw_bernoulli_iid(subvector, prob); }
}

void draw_bernoulli_iid(vector<indicator_t>& data, size_t count, double prob) {
    assert(count >= 0);
    assert(prob >= 0 and prob <= 1);
    if (data.size() > 0) { WARNING("Erasing non-empty vector"); }

    data = vector<indicator_t>(count, false);
    std::bernoulli_distribution distrib(prob);
    for (auto&& elem : data) { elem = distrib(Random::global_gen); }
}

void draw_bernoulli_iid(vector<vector<indicator_t>>& data, size_t count_x, size_t count_y,
    const vector<double>& probs) {
    /* -- */
    assert(count_x >= 0 and count_y >= 0);
    assert(count_x == probs.size());
    if (data.size() > 0) { WARNING("Erasing non-empty vector"); }

    data = vector<vector<indicator_t>>(count_x, vector<indicator_t>(count_y, false));
    for (size_t i = 0; i < data.size(); i++) { draw_bernoulli_iid(data.at(i), probs.at(i)); }
}

void draw_bernoulli_iid(
    vector<vector<indicator_t>>& data, size_t count_x, size_t count_y, double prob) {
    assert(count_x >= 0 and count_y >= 0);
    assert(prob >= 0 and prob <= 1);
    if (data.size() > 0) { WARNING("Erasing non-empty vector"); }

    data = vector<vector<indicator_t>>(count_x, vector<indicator_t>(count_y, false));
    for (auto&& subvector : data) { draw_bernoulli_iid(subvector, prob); }
}