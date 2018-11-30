#pragma once

#include "BufferManager.hpp"
#include "partition.hpp"

class PartitionedBufferManager {
    Partition partition;
    std::vector<BufferManager> managers;

  public:
    PartitionedBufferManager(Partition partition)
        : partition(partition), managers(partition.size()) {}

    template <class T>
    void add(T& x) {
        static_assert(is_partitionable<T>::value,
            "PartitionedBufferManager::add: type T is not partitionable");
        assert(x.size() % partition.size_all() == 0);

        size_t multiplicity = x.size() / partition.size_all();
        int i = 0;
        for (size_t subset = 0; subset < partition.size(); subset++) {
            size_t subset_size = partition.partition_size(subset);
            size_t nb_elements = subset_size * multiplicity;
            managers.at(subset).add_subset(x, i, nb_elements);
        }
    }
};