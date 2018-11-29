#pragma once

#include <mpi.h>
#include <cstdio>

class SendBuffer {
    size_t allocated_buffer_size{16};
    void* buffer{malloc(16)};

    size_t buffer_size{0}; // should be identical
    int buffer_position{0};

    size_t int_size() const {
        int result;
        MPI_Type_size(MPI_INT, &result);
        return static_cast<size_t>(result);
    }

    size_t double_size() const {
        int result;
        MPI_Type_size(MPI_DOUBLE, &result);
        return static_cast<size_t>(result);
    }

    void double_buffer_size() {
        void* new_buffer = malloc(allocated_buffer_size * 2);
        memcpy(new_buffer, buffer, buffer_size);
        free(buffer);
        allocated_buffer_size *= 2;
        buffer = new_buffer;
    }

  public:
    SendBuffer() {
        int packed_size;
        MPI_Type_size(MPI_PACKED, &packed_size);
        assert(packed_size == 1);
    }

    ~SendBuffer() { free(buffer); }

    void pack(int* data, size_t count = 1) {
        size_t add_size = int_size() * count;
        if (buffer_size + add_size > allocated_buffer_size) {
            double_buffer_size();
            pack(data, count);
        } else {
            assert(buffer_position == static_cast<int>(buffer_size));
            buffer_size += add_size;
            MPI_Pack(data, count, MPI_INT, buffer, allocated_buffer_size, &buffer_position, MPI_COMM_WORLD);
            assert(buffer_position == static_cast<int>(buffer_size));
        }
    }

    void pack(double* data, size_t count = 1) {
        size_t add_size = double_size() * count;
        if (buffer_size + add_size > allocated_buffer_size) {
            double_buffer_size();
            pack(data, count);
        } else {
            assert(buffer_position == static_cast<int>(buffer_size));
            buffer_size += add_size;
            MPI_Pack(data, count, MPI_DOUBLE, buffer, allocated_buffer_size, &buffer_position, MPI_COMM_WORLD);
            assert(buffer_position == static_cast<int>(buffer_size));
        }
    }
};

class Serializable {
  public:
    virtual size_t buffer_size() = 0;
    virtual void serialize(SendBuffer&) = 0;
};