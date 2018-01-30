
#ifndef MULTIPROBMODEL_H
#define MULTIPROBMODEL_H

#include "ProbModel.hpp"
#include "MultiGeneMPIModule.hpp"

class MultiGeneProbModel : public ProbModel, public MultiGeneMPIModule {

    public:

    MultiGeneProbModel(int inmyid, int innprocs) : ProbModel(), MultiGeneMPIModule(inmyid, innprocs) {}

    virtual void Update() override {
        if (! myid) {
            MasterUpdate();
        }
        else    {
            SlaveUpdate();
        }
    }

    virtual double Move() override {
        if (! myid) {
            MasterMove();
        }
        else    {
            SlaveMove();
        }
        return 1;
    }

	virtual void FromStream(istream& is) override {
        if (! myid) {
            MasterFromStream(is);
        }
        else    {
            SlaveFromStream();
        }
    }

    virtual void ToStream(ostream& os) const override {
        if (! myid) {
            MasterToStream(os);
        }
        else    {
            SlaveToStream();
        }
    }

    virtual void MasterToStream(ostream& os) const {}
    virtual void SlaveToStream() const {}
    virtual void MasterFromStream(istream& is) {}
    virtual void SlaveFromStream() {}

    virtual void MasterMove() {}
    virtual void SlaveMove() {}

    virtual void MasterUpdate() {}
    virtual void SlaveUpdate() {}
};

#endif

