#pragma once

#include <string>
#include "decl_info.hpp"
using std::string;

/*==================================================================================================
  Various utilities to help declare introspection interfaces
==================================================================================================*/

template <class User, class Target, class... Tags>
void declare(User& user, string name, Target& target) {
    user.template process_declaration<Tags...>(target);
}