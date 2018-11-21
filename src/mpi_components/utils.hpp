#pragma once

#include <assert.h>
#include <mpi.h>
#include <vector>

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

template <typename T>
MPI_Datatype get_datatype() {
    return GetMPIDatatype<T>()();
}

struct StructMetaData {
    std::vector<MPI_Datatype> types;
    std::vector<MPI_Aint> offsets;
    std::vector<int> block_lengths;

    MPI_Datatype declare_and_commit() {
        assert(types.size() == offsets.size());
        assert(types.size() == block_lengths.size());
        MPI_Datatype result;
        MPI_Type_create_struct(
            types.size(), block_lengths.data(), offsets.data(), types.data(), &result);
        MPI_Type_commit(&result);
        return result;
    }
};

#define STRUCT_DECL(CLASS)                          \
    StructMetaData struct_decl_##CLASS;             \
    StructMetaData& struct__ = struct_decl_##CLASS; \
    struct__.types.push_back(MPI_UB);               \
    struct__.offsets.push_back(sizeof(CLASS));      \
    struct__.block_lengths.push_back(0);            \
    CLASS x__;
#define ATTRIBUTE(NAME)                                           \
    struct__.types.push_back(get_datatype<decltype(x__.NAME)>()); \
    struct__.offsets.push_back(offsetof(decltype(x__), NAME));    \
    struct__.block_lengths.push_back(1);
#define STRUCT_DECL_COMMIT(NAME) NAME = struct__.declare_and_commit();