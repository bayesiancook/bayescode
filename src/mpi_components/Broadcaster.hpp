#pragma once

#include <functional>
#include "Process.hpp"

/*==================================================================================================
  Broadcaster
  An object responsible for broadcasting the values of specified fields from one process to others
==================================================================================================*/
template <typename T>
class BroadcasterMaster {
    Process& p;
    MPI_Datatype datatype;
    std::vector<T> buf;
    std::vector<std::function<void()>> writers;

  public:
    BroadcasterMaster(Process& p = *MPI::p) : p(p) {
        if (std::is_same<T, double>()) {
            datatype = MPI_DOUBLE;
        } else if (std::is_same<T, int>()) {
            datatype = MPI_INT;
        } else {
            p.message("Error: template parameter T is not a supported datatype!");
            exit(1);
        }
    }

    void add(T& target) {
        writers.push_back([&target, this]() { buf.push_back(target); });
    }

    void write_buffer() {
        buf.clear();
        for (auto writer : writers) { writer(); }
        for (auto e : buf) { std::cout << e << " "; }
        std::cout << "\n";
    }

    void broadcast() {
        write_buffer();
        MPI_Bcast(buf.data(), buf.size(), datatype, p.rank, MPI_COMM_WORLD);
    }
};

template <typename T>
class BroadcasterSlave {
    using buf_it = typename std::vector<T>::iterator;

    Process& p;
    MPI_Datatype datatype;
    std::vector<T> buf;
    std::vector<std::function<buf_it(buf_it)>> readers;

  public:
    BroadcasterSlave(Process& p = *MPI::p) : p(p) {
        if (std::is_same<T, double>()) {
            datatype = MPI_DOUBLE;
        } else if (std::is_same<T, int>()) {
            datatype = MPI_INT;
        } else {
            p.message("Error: template parameter T is not a supported datatype!");
            exit(1);
        }
    }

    void add(T& target) {
        readers.push_back([&target, this](buf_it it) {
            target = *it;
            return it + 1;
        });
        buf.push_back(-1);
    }

    void read_buffer() {
        buf_it it = buf.begin();
        for (auto reader : readers) { it = reader(it); }
    }

    void broadcast() {
        MPI_Bcast(buf.data(), buf.size(), datatype, 0, MPI_COMM_WORLD);
        read_buffer();
    }
};
