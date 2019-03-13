#pragma once

#include <vector>
#include "tags/decl_utils.hpp"

/*
====================================================================================================
  Serializable stuff
==================================================================================================*/
// clang-format off

// is supported by buffers directly (int and double)
template <class T> struct is_default_serializable:
    std::integral_constant<bool, std::is_same<T, int>::value or std::is_same<T, double>::value> {};

// implements template <T> void serialization_interface(T&)
template <class T> struct has_custom_serialization : std::false_type {};

// implements size_t size() or T& operator[]()
template <class T> struct has_size :            std::false_type {};
template <class T> struct has_access_operator : std::false_type {};
template <class T> struct value_type       { using value = void; };

// can be partitioned by buffer manager
template <class T> struct is_partitionable :
    std::integral_constant<bool, has_size<T>::value and has_access_operator<T>::value> {};

// can be serialized by contiguous chunks (eg, vector<double>)
template <class T> struct is_contiguously_serializable : std::integral_constant<bool,
    is_partitionable<T>::value and is_default_serializable<typename value_type<T>::value>::value
> {};

template <class T> struct is_nontrivial_vector : std::integral_constant<bool,
    is_partitionable<T>::value and not is_default_serializable<typename value_type<T>::value>::value
> {};

/*
====================================================================================================
  Vector traits
==================================================================================================*/
template <class T> struct has_size<std::vector<T>> :            std::true_type {};
template <class T> struct has_access_operator<std::vector<T>> : std::true_type {};
template <class T> struct value_type<std::vector<T>>         { using value = T; };
// clang-format on

template <class T>
struct external_interface<std::vector<T>> {
    template <class Info, class Target>
    static void declare_interface(Info info, Target& target) {
        for (size_t i = 0; i < target.size(); i++) { declare(info, std::to_string(i), target[i]); }
    }
};