#pragma once

#include <algorithm>
#include <cstdio>
#include "utils.hpp"

/*
====================================================================================================
  SendBuffer
  A buffer which can pack ints and doubles for sending over MPI
==================================================================================================*/
class SendBuffer {
    size_t allocated_buffer_size{0};
    void* buffer{nullptr};

    size_t buffer_size{0};  // should be identical to position
    int buffer_position{0};


    void double_buffer_size() {
        assert(buffer_size <= allocated_buffer_size);
        allocated_buffer_size = std::max(static_cast<size_t>(8), allocated_buffer_size * 2);
        buffer = realloc(buffer, allocated_buffer_size);
    }

  public:
    SendBuffer() { assert(MPI::packed_size() == 1); }
    SendBuffer(SendBuffer const&) = delete;
    SendBuffer& operator=(SendBuffer const&) = delete;
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