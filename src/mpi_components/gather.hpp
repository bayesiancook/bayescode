#pragma once

#include "interfaces.hpp"
#include "partition.hpp"

/*==================================================================================================
  GatherMaster
==================================================================================================*/
template <typename T>
class GatherMaster : public Proxy, public RegistrarBase<GatherMaster<T>> {
    using buf_it = typename std::vector<T>::iterator;

    Process& p;
    MPI_Datatype datatype;
    std::vector<T> buf;
    std::vector<int> displs, revcounts;
    std::vector<std::function<void(std::vector<buf_it>&)>> readers;

    friend RegistrarBase<GatherMaster<T>>;

    void register_element(std::string, std::vector<T>& target, const Partition& partition) {
        for (int i = 1; i < p.size; i++) {
            size_t p_size = partition.partition_size(i);
            revcounts.at(i) += p_size;
            for (int j = i + 1; j < p.size; j++) { displs.at(j) += p_size; }
        }
        for (int i = 0; i < partition.size_all(); i++) { buf.push_back(-1); }
        readers.push_back([&target, partition, this](std::vector<buf_it>& its) {
            for (int i = 1; i < p.size; i++) {
                buf_it end = it.at(i) + partition.partition_size(i);
                target = std::vector<T>(it.at(i), end);
                it.at(i) = end;
            }
        });
    }

    void read_buffer() {
        std::vector<buf_it> its{buf.begin()};
        for (int i = 1; i < p.size; i++) { its.push_back(buf.begin() + displs.at(i)); }
        for (auto reader : readers) { reader(its); }
    }

  public:
    GatherMaster(Process& p = *MPI::p)
        : p(p), datatype(get_datatype<T>()), displs(p.size, 0), revcounts(p.size, 0) {}

    void acquire() final {
        MPI_Gatherv(NULL, 0, datatype, buf.data(), revcounts.data(), displs.data(), datatype, 0,
            MPI_COMM_WORLD);
        read_buffer();
    }
};

/*==================================================================================================
  GatherSlave
==================================================================================================*/
template <typename T>
class GatherSlave : public Proxy, public RegistrarBase<GatherSlave<T>> {
    Process& p;
    MPI_Datatype datatype;
    std::vector<T> buf;
    std::vector<std::function<void()>> writers;

    friend RegistrarBase<GatherSlave<T>>;

    void register_element(std::string, std::vector<T>& target) {
        writers.push_back(
            [&target, this]() { buf.insert(buf.end(), target.begin(), target.end()); });
    }

    void write_buffer() {
        buf.clear();
        for (auto writer : writers) { writer(); }
    }

  public:
    GatherSlave(Process& p = *MPI::p)
        : p(p), datatype(get_datatype<T>()), displs(p.size, 0), revcounts(p.size, 0) {}

    void release() final {
        write_buffer();
        MPI_Gatherv(
            buf.data(), buf.size(), datatype, NULL, NULL, NULL, datatype, 0, MPI_COMM_WORLD);
    }
};
