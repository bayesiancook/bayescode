#pragma once

#include <numeric>
#include "ReceiveBuffer.hpp"
#include "SendBuffer.hpp"

/*
====================================================================================================
  BufferManager
  A class to which ints and doubles (and vectors of ints and doubles) can be registered.
  Registered objects can then be packed in bulk into a buffer or read in bulk from a buffer.
==================================================================================================*/
class BufferManager {
    // clang-format off
    struct double_array_t { double* data; size_t size; };
    struct int_array_t    { int*    data; size_t size; };
    // clang-format on

    std::vector<int*> ints;
    std::vector<double*> doubles;
    std::vector<int_array_t> int_arrays;
    std::vector<double_array_t> double_arrays;

    std::unique_ptr<SendBuffer> _send_buffer;
    std::unique_ptr<ReceiveBuffer> _receive_buffer;

  public:
    BufferManager() = default;

    void add(int& x) { ints.push_back(&x); }
    void add(double& x) { doubles.push_back(&x); }
    void add(std::vector<int>& x) { int_arrays.push_back({x.data(), x.size()}); }
    void add(std::vector<double>& x) { double_arrays.push_back({x.data(), x.size()}); }

    template <class T>
    void add(T& x) {
        x.template serialization_interface<BufferManager>(*this);
    }

    template <class Arg, class... Args>
    void add(Arg& arg, Args&&... args) {
        add(arg);
        add(std::forward<Args>(args)...);
    }

    size_t nb_ints() const {
        auto sum_size = [](int acc, int_array_t v) { return acc + v.size; };
        int vec_size = std::accumulate(int_arrays.begin(), int_arrays.end(), 0, sum_size);
        return ints.size() + vec_size;
    }

    size_t nb_doubles() const {
        auto sum_size = [](int acc, double_array_t v) { return acc + v.size; };
        int vec_size = std::accumulate(double_arrays.begin(), double_arrays.end(), 0, sum_size);
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
        for (auto x : int_arrays) { _send_buffer->pack(x.data, x.size); }
        for (auto x : double_arrays) { _send_buffer->pack(x.data, x.size); }
        assert(buffer_size() == _send_buffer->size());
        return _send_buffer->data();
    }

    void receive() {
        assert(_receive_buffer.get() != nullptr);
        assert(_receive_buffer->size() == buffer_size());
        for (auto x : ints) { *x = _receive_buffer->unpack<int>(); }
        for (auto x : doubles) { *x = _receive_buffer->unpack<double>(); }
        for (auto x : int_arrays) { _receive_buffer->unpack_array<int>(x.data, x.size); }
        for (auto x : double_arrays) { _receive_buffer->unpack_array<double>(x.data, x.size); }
    }
};