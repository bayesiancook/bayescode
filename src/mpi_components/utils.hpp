#pragma once

#include <assert.h>
#include <mpi.h>
#include <vector>
#include "Process.hpp"
#include "components/RegistrarBase.hpp"
#include "global/logging.hpp"
#include "operations/interfaces.hpp"

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
        fprintf(stderr, "Error in GetMPIDatatype: type %s is not a registered MPI datatype!\n",
            demangle(typeid(T).name()).c_str());
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
  Some helpers to get the size of MPI types */

namespace MPI {
    size_t int_size() {
        int result;
        MPI_Type_size(MPI_INT, &result);
        return static_cast<size_t>(result);
    }

    size_t double_size() {
        int result;
        MPI_Type_size(MPI_DOUBLE, &result);
        return static_cast<size_t>(result);
    }

    size_t packed_size() {
        int result;
        MPI_Type_size(MPI_PACKED, &result);
        return static_cast<size_t>(result);
    }
};  // namespace MPI

/*--------------------------------------------------------------------------------------------------
  Function to get type with better syntax */

template <typename T>
MPI_Datatype get_datatype() {
    return GetMPIDatatype<T>()();
}

/*
====================================================================================================
  Group communication instantiation helpers
  Functions to help declare things globally and instantiate the correct objects
==================================================================================================*/

template <class Model, class Class>
using decl_pointer_t = void (Model::*)(RegistrarBase<Class>&);

template <class MasterClass, class SlaveClass, class Model, class... ConstructorArgs>
std::unique_ptr<Proxy> instantiate_and_declare(Model& m,
    decl_pointer_t<Model, MasterClass> f_master, decl_pointer_t<Model, SlaveClass> f_slave,
    filter_t filter, ConstructorArgs&&... args) {
    /* -- */
    std::unique_ptr<Proxy> result{nullptr};
    if (!MPI::p->rank) {
        auto component = new MasterClass(std::forward<ConstructorArgs>(args)...);
        component->register_from_method(m, f_master, filter);
        result.reset(dynamic_cast<Proxy*>(component));
    } else {
        auto component = new SlaveClass(std::forward<ConstructorArgs>(args)...);
        component->register_from_method(m, f_slave, filter);
        result.reset(dynamic_cast<Proxy*>(component));
    }
    return result;
}

template <class MasterClass, class SlaveClass, class Model, class... ConstructorArgs>
std::unique_ptr<Proxy> instantiate_and_declare_from_model(
    Model& m, filter_t filter, ConstructorArgs&&... args) {
    /* -- */
    return instantiate_and_declare<MasterClass, SlaveClass, Model>(m, &Model::declare_model,
        &Model::declare_model, filter, std::forward<ConstructorArgs>(args)...);
}

/*
====================================================================================================
  mpi_run
  Wrapper around a main-like function that initializes MPI and sets MPI::p to correspond to the
  local MPI process
==================================================================================================*/
template <class F, class... Args>
void mpi_run(int argc, char** argv, F f) {
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI::p = std::unique_ptr<Process>(new Process(rank, size));
    MPI::p->message("Started MPI process");
    f(argc, argv);
    MPI::p->message("End of MPI process");
    MPI_Finalize();
}
