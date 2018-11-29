#pragma once

#include <cstdio>
#include "utils.hpp"

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

class SendBuffer {
    size_t allocated_buffer_size{16};
    void* buffer{malloc(16)};

    size_t buffer_size{0};  // should be identical to position
    int buffer_position{0};


    void double_buffer_size() {
        void* new_buffer = malloc(allocated_buffer_size * 2);
        memcpy(new_buffer, buffer, buffer_size);
        free(buffer);
        allocated_buffer_size *= 2;
        buffer = new_buffer;
    }

  public:
    SendBuffer() { assert(MPI::packed_size() == 1); }

    ~SendBuffer() { free(buffer); }

    void pack(int* data, size_t count = 1) {
        size_t add_size = MPI::int_size() * count;
        if (buffer_size + add_size > allocated_buffer_size) {
            double_buffer_size();
            pack(data, count);
        } else {
            assert(buffer_position == static_cast<int>(buffer_size));
            buffer_size += add_size;
            MPI_Pack(data, count, MPI_INT, buffer, allocated_buffer_size, &buffer_position,
                MPI_COMM_WORLD);
            assert(buffer_position == static_cast<int>(buffer_size));
        }
    }

    void pack(double* data, size_t count = 1) {
        size_t add_size = MPI::double_size() * count;
        if (buffer_size + add_size > allocated_buffer_size) {
            double_buffer_size();
            pack(data, count);
        } else {
            assert(buffer_position == static_cast<int>(buffer_size));
            buffer_size += add_size;
            MPI_Pack(data, count, MPI_DOUBLE, buffer, allocated_buffer_size, &buffer_position,
                MPI_COMM_WORLD);
            assert(buffer_position == static_cast<int>(buffer_size));
        }
    }

    size_t size() const { return buffer_size; }

    void* data() const { return buffer; }
};

class ReceiveBuffer {
    void* buffer{nullptr};
    int buffer_position{0};
    size_t buffer_size{0};

  public:
    ReceiveBuffer(void* buffer, size_t buffer_size) : buffer(buffer), buffer_size(buffer_size) {}

    template <typename T>
    T unpack() {
        T result;
        MPI_Datatype mpi_type = GetMPIDatatype<T>()();
        MPI_Unpack(buffer, buffer_size, &buffer_position, &result, 1, mpi_type, MPI_COMM_WORLD);
        return result;
    }
};