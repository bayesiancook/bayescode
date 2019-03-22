#pragma once

#include <assert.h>
#include <vector>
#include "Random.hpp"
using std::vector;

/*--------------------------------------------------------------------------------------------------
  Type aliases to make delcarations more informative */
template <class T>
using cond_vector = vector<T>;

template <class T>
using site_vector = vector<T>;

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
void set_all_to(vector<T>& data, int count, T value) {
    assert(count >= 0);
    if (data.size() > 0) { WARNING("Erasing non-empty vector"); }
    data = vector<T>(count, value);
}

template <class T>
void set_all_to(vector<vector<T>>& data, int count_x, int count_y, T value) {
    assert(count >= 0);
    if (data.size() > 0) { WARNING("Erasing non-empty vector"); }
    data = vector<vector<T>>(count_x, vector<T>(count_y, value));
}

/*--------------------------------------------------------------------------------------------------
  Counting bools set to true in vectors */
int count_indicators(vector<bool>& data) {
    int result = 0;
    for (auto ind : data) { result += ind ? 1 : 0; }
    return result;
}

int count_indicators(vector<vector<bool>>& data) {
    int result = 0;
    for (auto&& subvector : data) { result += count_indicators(subvector); }
    return result;
}

/*--------------------------------------------------------------------------------------------------
  Drawing vector of indicators from iid bernoullis */
void draw_bernoulli_iid(vector<bool>& data, double prob) {
    assert(prob >= 0 and prob <= 1);
    std::bernoulli_distribution distrib(prob);
    for (auto&& elem : data) { elem = distrib(Random::global_gen); }
}

void draw_bernoulli_iid(vector<vector<bool>>& data, const vector<double>& probs) {
    assert(probs.size() == data.size());
    for (size_t i = 0; i < data.size(); i++) { draw_bernoulli_iid(data.at(i), probs.at(i)); }
}

void draw_bernoulli_iid(vector<vector<bool>>& data, double prob) {
    assert(prob >= 0 and prob <= 1);
    for (auto&& subvector : data) { draw_bernoulli_iid(subvector, prob); }
}

void draw_bernoulli_iid(vector<bool>& data, int count, double prob) {
    assert(count >= 0);
    assert(prob >= 0 and prob <= 1);
    if (data.size() > 0) { WARNING("Erasing non-empty vector"); }
    
    data = vector<bool>(count, false);
    std::bernoulli_distribution distrib(prob);
    for (auto&& elem : data) { elem = distrib(Random::global_gen); }
}

void draw_bernoulli_iid(
    vector<vector<bool>>& data, int count_x, int count_y, const vector<double>& probs) {
    /* -- */
    assert(count_x >= 0 and count_y >= 0);
    assert(probs.size() == data.size());
    if (data.size() > 0) { WARNING("Erasing non-empty vector"); }
    
    data = vector<vector<bool>>(count_x, vector<bool>(count_y, false));
    for (size_t i = 0; i < data.size(); i++) { draw_bernoulli_iid(data.at(i), probs.at(i)); }
}

void draw_bernoulli_iid(vector<vector<bool>>& data, int count_x, int count_y, double prob) {
    assert(count_x >= 0 and count_y >= 0);
    assert(prob >= 0 and prob <= 1);
    if (data.size() > 0) { WARNING("Erasing non-empty vector"); }
    
    data = vector<vector<bool>>(count_x, vector<bool>(count_y, false));
    for (auto&& subvector : data) { draw_bernoulli_iid(subvector, prob); }
}