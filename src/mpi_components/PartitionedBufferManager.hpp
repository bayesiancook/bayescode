#pragma once

#include "BufferManager.hpp"
#include "partition.hpp"

class PartitionedBufferManager {
    Partition partition;
    std::vector<BufferManager> temp_managers;

    BufferManager manager;
    bool manager_ready{false};

    void check_manager() {
        if (!manager_ready) {
            for (size_t subset = 0; subset < partition.size(); subset++) {
                manager.merge(temp_managers.at(subset));
            }
        }
    }

  public:
    PartitionedBufferManager(Partition partition)
        : partition(partition), temp_managers(partition.size()) {}

    template <class T>
    void add(T& x) {
        static_assert(is_partitionable<T>::value,
            "PartitionedBufferManager::add: type T is not partitionable");
        assert(x.size() % partition.size_all() == 0);
        assert(!manager_ready);

        size_t multiplicity = x.size() / partition.size_all();
        int i = 0;
        for (size_t subset = 0; subset < partition.size(); subset++) {
            size_t subset_size = partition.partition_size(subset);
            size_t nb_elements = subset_size * multiplicity;
            temp_managers.at(subset).add_subset(x, i, nb_elements);
            i += nb_elements;
        }
    }

    void* send_buffer() {
        check_manager();
        return manager.send_buffer();
    }

    void* receive_buffer() {
        check_manager();
        return manager.receive_buffer();
    }

    void receive() {
        check_manager();
        manager.receive();
    }

    size_t buffer_size() {
        check_manager();
        return manager.buffer_size();
    }
};