#pragma once

#include "ParamConfig.hpp"

template <class ValueT, class HyperValueT>
istream &operator>>(istream &is, param::Config<ValueT, HyperValueT> &c) {
    auto contains = [](string input, char sep) {
        return std::find(input.begin(), input.end(), sep) != input.end();
    };

    auto split_at = [](string input, char sep) {
        auto split_point = std::find(input.begin(), input.end(), sep);
        assert(split_point != input.end());
        return make_pair(string(input.begin(), split_point), string(split_point + 1, input.end()));
    };

    string declaration;
    is >> declaration;

    if (contains(declaration, ':')) {
        assert(not contains(declaration, '/'));
        auto colon_split = split_at(declaration, ':');
        std::stringstream(colon_split.first) >> c.mode;
        assert(c.mode == param::shared or c.mode == param::fixed);
        std::stringstream(colon_split.second) >> c.value;
    } else if (contains(declaration, '/')) {
        auto slash_split = split_at(declaration, '/');
        std::stringstream(slash_split.first) >> c.mode;
        assert(c.mode == param::independent or c.mode == param::shrunk);
        std::stringstream(slash_split.second) >> c.hyper;
    } else {
        std::stringstream(declaration) >> c.mode;
        assert(c.mode == param::shared or c.mode == param::shrunk);
    }
    return is;
}