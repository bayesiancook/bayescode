#pragma once

#include <vector>

/*
====================================================================================================
  Serializable stuff
==================================================================================================*/
// clang-format off

// is supported by buffers directly (int and double)
template <class T> struct is_default_serializable:
    std::integral_constant<bool, std::is_same<T, int>::value or std::is_same<T, double>::value> {};

// implements template <T> void serialization_interface(T&)
template <class T> struct has_custom_serialization :     public std::false_type {};

// can be serialized by contiguous chunks (eg, vector<double>)
template <class T> struct is_contiguously_serializable : public std::false_type {};

/*
====================================================================================================
  Array stuff
==================================================================================================*/

// implements size_t size() or T& operator[]()
template <class T> struct has_size :            public std::false_type {};
template <class T> struct has_access_operator : public std::false_type {};

// can be partitioned by buffer manager
template <class T> struct is_partitionable :
    std::integral_constant<bool, has_size<T>::value and has_access_operator<T>::value> {};

/*
====================================================================================================
  Vector traits
==================================================================================================*/
template <class T> struct is_contiguously_serializable<std::vector<T>>
    { static bool const value = is_default_serializable<T>::value; };

template <class T> struct has_size<std::vector<T>> :            public std::true_type {};
template <class T> struct has_access_operator<std::vector<T>> : public std::true_type {};
// clang-format on