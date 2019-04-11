#pragma once

#include "global/Random.hpp"
#include "global/defs.hpp"
#include "type_aliases.hpp"

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