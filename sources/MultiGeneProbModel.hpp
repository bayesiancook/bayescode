
#ifndef MULTIPROBMODEL_H
#define MULTIPROBMODEL_H

#include "ProbModel.hpp"
#include "MultiGeneMPIModule.hpp"

class MultiGeneProbModel : public ProbModel, public MultiGeneMPIModule {

    public:

    MultiGeneProbModel(int inmyid, int innprocs) : ProbModel(), MultiGeneMPIModule(inmyid, innprocs) {}

    virtual double Move() override {
        if (! myid) {
            MasterMove();
        }
        else    {
            SlaveMove();
        }
        return 1;
    }

    virtual void MasterMove() {}
    virtual void SlaveMove() {}

};

#endif

