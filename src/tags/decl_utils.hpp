#pragma once

#include <string>
#include "decl_info.hpp"
#include "processing_bricks.hpp"
using std::string;

/*==================================================================================================
  Various utilities to help declare introspection interfaces and apply them
==================================================================================================*/

/*--------------------------------------------------------------------------------------------------
  Basic declare function to be used in declare_interface methods */
template <class... Tags, class PrInfo, class... Args>
void declare(PrInfo& processing_info, string name, Args&&... args) {
    auto decl_info = make_decl_info<Tags...>(name);
    PrInfo::processing::forward_declaration(
        processing_info, decl_info, std::forward<Args>(args)...);
}

/*--------------------------------------------------------------------------------------------------
  Basic apply: apply interface declaration of provider to user */
template <class User, class Provider>
void basic_apply(User& user, Provider& provider) {
    auto processing_info = make_processing_info<processing::End>(user);
    provider.declare_interface(processing_info);
}

/*--------------------------------------------------------------------------------------------------
  Filter apply: allows application of only declarations with a given tag */
template <class Tag, class User, class Provider>
void filter_apply(User& user, Provider& provider) {
    using namespace processing;
    auto processing_info = make_processing_info<FilterTag<Tag, End>>(user);
    provider.declare_interface(processing_info);
}

/*--------------------------------------------------------------------------------------------------
  Typefilter apply: allows application of only declarations with a given target type */
template <class Type, class User, class Provider>
void typefilter_apply(User& user, Provider& provider) {
    using namespace processing;
    auto processing_info = make_processing_info<FilterType<Type, End>>(user);
    provider.declare_interface(processing_info);
}

/*--------------------------------------------------------------------------------------------------
  Single unroll of structures with tag Tag */
template <class Tag, class User, class Provider>
void unrollif_apply(User& user, Provider& provider) {
    using namespace processing;
    auto processing_info = make_processing_info<UnrollIf<Tag, End>>(user);
    provider.declare_interface(processing_info);
}

/*--------------------------------------------------------------------------------------------------
  Recursive unroll of structures with tag Tag */
template <class Tag, class User, class Provider>
void recif_apply(User& user, Provider& provider) {
    using namespace processing;
    auto processing_info = make_processing_info<RecursiveUnroll<HasTag<Tag>, End>>(user);
    provider.declare_interface(processing_info);
}