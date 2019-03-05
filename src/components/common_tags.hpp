// TODO rename file

#pragma once
#include "tags/decl_utils.hpp"

// to be used in declarations to indicate a given declaration is a sub structure
struct SubStructure {};

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