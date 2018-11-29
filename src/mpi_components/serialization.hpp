#pragma once

#include <algorithm>
#include <cstdio>
#include "utils.hpp"

/*
====================================================================================================
  Some helpers to get the size of MPI types
==================================================================================================*/
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
};

/*
====================================================================================================
  BufferManager
  A class to which ints and doubles (and vectors of ints and doubles) can be registered.
  Registered objects can then be packed in bulk into a buffer or read in bulk from a buffer.
==================================================================================================*/
class BufferManager {
    std::vector<int*> ints;
    std::vector<double*> doubles;
    std::vector<std::vector<int>*> int_vectors;
    std::vector<std::vector<double>*> double_vectors;

    std::unique_ptr<SendBuffer> _send_buffer;
    std::unique_ptr<ReceiveBuffer> _receive_buffer;

  public:
    BufferManager() = default;

    void add(int& x) { ints.push_back(&x); }
    void add(double& x) { doubles.push_back(&x); }
    void add(std::vector<int>& x) { int_vectors.push_back(&x); }
    void add(std::vector<double>& x) { double_vectors.push_back(&x); }

    template <class T>
    void add(T& x) {
        x.template serialization_interface<BufferManager>(*this);
    }

    size_t nb_ints() const {
        auto sum_size = [](int acc, std::vector<int>* v) { return acc + v->size(); };
        int vec_size = std::accumulate(int_vectors.begin(), int_vectors.end(), 0, sum_size);
        return ints.size() + vec_size;
    }

    size_t nb_doubles() const {
        auto sum_size = [](int acc, std::vector<double>* v) { return acc + v->size(); };
        int vec_size = std::accumulate(double_vectors.begin(), double_vectors.end(), 0, sum_size);
        return doubles.size() + vec_size;
    }

    size_t buffer_size() const {
        return nb_ints() * MPI::int_size() + nb_doubles() * MPI::double_size();
    }

    void* receive_buffer() {
        _receive_buffer.reset(new ReceiveBuffer(buffer_size()));
        return _receive_buffer->data();
    }

    void* send_buffer() {
        _send_buffer.reset(new SendBuffer());
        for (auto x : ints) { _send_buffer->pack(x); }
        for (auto x : doubles) { _send_buffer->pack(x); }
        for (auto x : int_vectors) { _send_buffer->pack(x->data(), x->size()); }
        for (auto x : double_vectors) { _send_buffer->pack(x->data(), x->size()); }
        assert(buffer_size() == _send_buffer->size());
        return _send_buffer->data();
    }

    // second param is optional but can be used to ensure size is correct
    void receive() {
        assert(_receive_buffer.get() != nullptr);
        assert(_receive_buffer->size() == buffer_size());
        for (auto x : ints) { *x = _receive_buffer->unpack<int>(); }
        for (auto x : doubles) { *x = _receive_buffer->unpack<double>(); }
        for (auto x : int_vectors) { *x = _receive_buffer->unpack_vector<int>(x->size()); }
        for (auto x : double_vectors) { *x = _receive_buffer->unpack_vector<double>(x->size()); }
    }
};