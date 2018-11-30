#pragma once

#include <cstdio>
#include <vector>
#include "utils.hpp"

/*
====================================================================================================
  ReceiveBuffer
  A class to unpack ints and doubles from a MPI buffer
==================================================================================================*/
class ReceiveBuffer {
    void* buffer{nullptr};
    int buffer_position{0};
    size_t buffer_size{0};

  public:
    ReceiveBuffer(size_t buffer_size = 0) : buffer_size(buffer_size) {
        buffer = malloc(buffer_size);
    }
    ReceiveBuffer(ReceiveBuffer const&) = delete;
    ReceiveBuffer& operator=(ReceiveBuffer const&) = delete;
    ~ReceiveBuffer() { free(buffer); }

    void* data() const { return buffer; }

    size_t size() const { return buffer_size; }

    template <typename T>
    T unpack() {
        T result;
        MPI_Datatype mpi_type = get_datatype<T>();
        MPI_Unpack(buffer, buffer_size, &buffer_position, &result, 1, mpi_type, MPI_COMM_WORLD);
        assert(buffer_position <= static_cast<int>(buffer_size));
        return result;
    }

    template <typename T>
    std::vector<T> unpack_vector(int count) {
        std::vector<T> result(count);
        MPI_Datatype mpi_type = get_datatype<T>();
        MPI_Unpack(
            buffer, buffer_size, &buffer_position, result.data(), count, mpi_type, MPI_COMM_WORLD);
        assert(buffer_position <= static_cast<int>(buffer_size));
        return result;
    }

    template <typename T>
    void unpack_array(T* data, size_t size) {
        MPI_Datatype mpi_type = get_datatype<T>();
        MPI_Unpack(buffer, buffer_size, &buffer_position, data, size, mpi_type, MPI_COMM_WORLD);
        assert(buffer_position <= static_cast<int>(buffer_size));
    }
};