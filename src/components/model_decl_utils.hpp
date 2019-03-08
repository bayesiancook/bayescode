#pragma once
#include "tags/decl_utils.hpp"

struct ModelNode {};

struct Stat {};

template <class... Tags, class... Args>
void model_node(Args&&... args) {
    declare<ModelNode, Tags...>(std::forward<Args>(args)...);
}

template <class... Tags, class... Args>
void model_stat(Args&&... args) {
    declare<Stat, Tags...>(std::forward<Args>(args)...);
}