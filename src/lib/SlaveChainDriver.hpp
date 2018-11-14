#pragma once
#include <string>
#include <vector>

class SlaveChainDriver {
    static std::string get_first_token(std::istream& is) {
        std::string result;
        is >> result;
        return result;
    }

    static int SlaveReceiveRunningStatus() {
        int status;
        MPI_Bcast(&status, 1, MPI_INT, 0, MPI_COMM_WORLD);
        return status;
    }

  public:
    SlaveChainDriver(std::string name, int every, int until, int size = 0)
        : name(name), every(every), until(until), size(size) {}

    void go() {
        if (size == 0)
            for (auto c : components) c->start();
        while (SlaveReceiveRunningStatus()) {
            for (int i = 0; i < every; i++)
                for (auto c : components) c->move(size * every + i);
            for (auto c : components) c->savepoint(size);
            size++;
        }
        for (auto c : components) c->end();
    }

    void add(ChainComponent& component) { components.push_back(&component); }

    SlaveChainDriver(std::istream& is) : name(get_first_token(is)) {
        is >> every;
        is >> until;
        is >> size;
    }

  private:
    std::string name;
    std::vector<ChainComponent*> components;
    //! saving frequency (i.e. number of move cycles performed between each point
    //! saved to file)
    int every{1};
    //! intended final size of the chain (until==-1 means no a priori specified
    //! upper limit)
    int until{-1};
    //! current size (number of points saved to file)
    int size{0};
};
