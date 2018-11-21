#pragma once

#include <assert.h>
#include <mpi.h>
#include <vector>
#include "interfaces.hpp"
template <class T>
class RegistrarBase;

/*
====================================================================================================
  Type Helpers
  Various things to help get MPI types from regular types
==================================================================================================*/

/*--------------------------------------------------------------------------------------------------
  GetMPIDatatype overloads */

template <class T>
struct GetMPIDatatype {
    MPI_Datatype operator()() {
        fprintf(
            stderr, "Error in get_datatype: template parameter T is not a supported datatype!\n");
        exit(1);
    }
};

template <>
struct GetMPIDatatype<double> {
    MPI_Datatype operator()() { return MPI_DOUBLE; }
};

template <>
struct GetMPIDatatype<float> {
    MPI_Datatype operator()() { return MPI_FLOAT; }
};

template <>
struct GetMPIDatatype<int> {
    MPI_Datatype operator()() { return MPI_INT; }
};

/*--------------------------------------------------------------------------------------------------
  Function to get type with better syntax */

template <typename T>
MPI_Datatype get_datatype() {
    return GetMPIDatatype<T>()();
}

struct StructMetaData {
    std::vector<MPI_Datatype> types;
    std::vector<MPI_Aint> offsets;
    std::vector<int> block_lengths;
    MPI_Aint size;

    MPI_Datatype declare_and_commit() {
        assert(types.size() == offsets.size());
        assert(types.size() == block_lengths.size());
        MPI_Datatype tmp;
        MPI_Type_create_struct(
            types.size(), block_lengths.data(), offsets.data(), types.data(), &tmp);
        MPI_Datatype result;
        MPI_Type_create_resized(tmp, 0, size, &result);
        MPI_Type_commit(&result);
        return result;
    }
};

/*--------------------------------------------------------------------------------------------------
  Macros to declare mpi sruct types */

#define STRUCT_DECL(CLASS)                          \
    StructMetaData struct_decl_##CLASS;             \
    StructMetaData& struct__ = struct_decl_##CLASS; \
    struct__.size = sizeof(CLASS);                  \
    CLASS x__;
#define ATTRIBUTE(NAME)                                           \
    struct__.types.push_back(get_datatype<decltype(x__.NAME)>()); \
    struct__.offsets.push_back(offsetof(decltype(x__), NAME));    \
    struct__.block_lengths.push_back(1);
#define STRUCT_COMMIT(NAME) NAME = struct__.declare_and_commit();
#define STRUCT_GLOBAL_DECL(CLASS, TYPENAME)            \
    MPI_Datatype TYPENAME;                             \
    template <>                                        \
    struct GetMPIDatatype<CLASS> {                     \
        MPI_Datatype operator()() { return TYPENAME; } \
    };

/*
====================================================================================================
  Group communication instantiation helpers
  Functions to help declare things globally and instantiate the correct objects
==================================================================================================*/

template <class Model, class Class>
using decl_pointer_t = void (Model::*)(RegistrarBase<Class>&);

using filter_t = std::set<std::string>;

template <class MasterClass, class SlaveClass, class Model>
std::unique_ptr<Proxy> instantiate_and_declare(Model& m,
    decl_pointer_t<Model, MasterClass> f_master, decl_pointer_t<Model, SlaveClass> f_slave,
    filter_t filter) {
    /* -- */
    assert(filter.size() != 0);
    std::unique_ptr<Proxy> result{nullptr};
    if (!MPI::p->rank) {
        auto component = new MasterClass();
        component->register_from_method(m, f_master, filter);
        result.reset(dynamic_cast<Proxy*>(component));
    } else {
        auto component = new SlaveClass();
        component->register_from_method(m, f_slave, filter);
        result.reset(dynamic_cast<Proxy*>(component));
    }
    return result;
}

template <class MasterClass, class SlaveClass, class Model>
std::unique_ptr<Proxy> instantiate_and_declare_from_model(Model& m, filter_t filter) {
    /* -- */
    return instantiate_and_declare<MasterClass, SlaveClass, Model>(
        m, &Model::declare_model, &Model::declare_model, filter);
}
