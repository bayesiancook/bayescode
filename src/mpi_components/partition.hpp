#pragma once

#include <map>
#include <numeric>
#include <set>
#include <string>
#include <vector>
#include "Process.hpp"

/*==================================================================================================
  Data types
==================================================================================================*/
using Index = std::string;
using IndexSet = std::set<Index>;
using IndexMapping = std::map<Index, Index>;

/*==================================================================================================
  Partition
==================================================================================================*/
class Partition {
    std::map<int, IndexSet> partition;
    const Process& process;  // local process (used in my* member functions)

  public:
    Partition(IndexSet indexes, size_t size, size_t first = 0, const Process& process = *MPI::p)
        : process(process) {
        size_t nb_indexes = indexes.size();
        for (size_t i = 0; i < size; i++) {
            auto begin = indexes.begin();
            std::advance(begin, i * nb_indexes / size);
            auto end = indexes.begin();
            std::advance(end, (i + 1) * nb_indexes / size);
            partition.insert({i + first, IndexSet(begin, end)});
        }
    }

    /*==============================================================================================
      Partition getters */
    IndexSet get_partition(int i) const {
        if (partition.find(i) != partition.end()) {
            return partition.at(i);
        } else {
            MPI::p->message("Error in get_partition: no subset for index %d", i);
            exit(1);
        }
    }

    IndexSet get_all() const {
        IndexSet result;
        for (auto subpartition : partition) {
            result.insert(subpartition.second.begin(), subpartition.second.end());
        }
        return result;
    }

    IndexSet my_partition() const { return get_partition(process.rank); }

    /*==============================================================================================
      Size information */
    size_t partition_size(int i) const {
        if (partition.find(i) != partition.end()) {
            return partition.at(i).size();
        } else {
            return size_all();
            // MPI::p->message("Error in partition_size: no subset for index %d", i);
            // exit(1);
        }
    }

    size_t my_partition_size() const { return partition_size(process.rank); }

    size_t size_all() const {
        size_t result = 0;
        for (auto subset : partition) { result += subset.second.size(); }
        return result;
    }

    size_t size() const { return partition.size(); }

    size_t max_partition_size() const {
        size_t result = 0;
        for (auto subset : partition) {
            result = subset.second.size() > result ? subset.second.size() : result;
        }
        return result;
    }

    /*==============================================================================================
      Other */
    int owner(Index index) const {
        for (auto subset : partition) {
            if (subset.second.find(index) != subset.second.end()) { return subset.first; }
        }
        return -1;
    }
};