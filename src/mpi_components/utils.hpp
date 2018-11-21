#pragma once

#include <mpi.h>
#include <vector>

#define STRUCT_DECL(CLASS)                          \
    StructMetaData struct_decl_##CLASS;             \
    StructMetaData& struct__ = struct_decl_##CLASS; \
    struct__.member_types.push_back(MPI_UB);        \
    struct__.offsets.push_back(sizeof(CLASS));      \
    CLASS x__;
#define ATTRIBUTE(name)                                                  \
    struct__.member_types.push_back(get_datatype<decltype(x__.name)>()); \
    struct__.offsets.push_back(offsetof(decltype(x__), name));


struct StructMetaData {
    std::vector<MPI_Datatype> member_types;
    std::vector<MPI_Aint> offsets;
};

template <typename T>
MPI_Datatype get_datatype() {
    if (std::is_same<T, double>()) {
        return MPI_DOUBLE;
    } else if (std::is_same<T, int>()) {
        return MPI_INT;
    } else {
        fprintf(
            stderr, "Error in get_datatype: template parameter T is not a supported datatype!\n");
        exit(1);
    }
}