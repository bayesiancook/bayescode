#ifndef MULTIGENE_H
#define MULTIGENE_H

#include <vector>
#include <string>
using namespace std;
#include "SequenceAlignment.hpp"
#include "Array.hpp"
#include "BranchArray.hpp"
#include "MPIBuffer.hpp"
#include "Parallel.hpp"

class MultiGeneMPIModule    {

    public:

    MultiGeneMPIModule(int inmyid, int innprocs) : myid(inmyid), nprocs(innprocs) {}
    ~MultiGeneMPIModule() {}

	int GetMyid() const {
		return myid;
	}

	int GetNprocs() const {
		return nprocs;
	}

    int GetNgene() const {
        return Ngene;
    }

    int GetLocalNgene() const   {
        return LocalNgene;
    }

    int GetLocalTotNsite() const    {
        return SlaveTotNsite[myid];
    }

    int GetTotNsite() const {
        return SlaveTotNsite[0];
    }

    int GetSlaveNgene(int proc) const   {
        if (myid)   {
            cerr << "error: slave in GetSlaveNgene\n";
            exit(1);
        }
        return SlaveNgene[proc];
    }

    string GetLocalGeneName(int gene) const {
        return GeneName[gene];
    }

    int GetLocalGeneNsite(int gene) const   {
        return GeneNsite[gene];
    }

    int GetSlaveTotNsite(int proc) const    {
        return SlaveTotNsite[proc];
    }

    void AllocateAlignments(string datafile);

    void PrintGeneList(ostream& os) const;

    template<class T> void MasterSendGlobal(const T& t) const   {
        MPIBuffer buffer(MPISize(t));
        buffer << t;
        MPI_Bcast(buffer.GetBuffer(),buffer.GetSize(),MPI_DOUBLE,0,MPI_COMM_WORLD);
    }

    template<class T> void SlaveReceiveGlobal(T& t) {
        MPIBuffer buffer(MPISize(t));
        MPI_Bcast(buffer.GetBuffer(),buffer.GetSize(),MPI_DOUBLE,0,MPI_COMM_WORLD);
        buffer >> t;
    }

    template<class T, class U> void MasterSendGlobal(const T& t, const U& u) const  {
        MPIBuffer buffer(MPISize(t) + MPISize(u));
        buffer << t << u;
        MPI_Bcast(buffer.GetBuffer(),buffer.GetSize(),MPI_DOUBLE,0,MPI_COMM_WORLD);
    }

    template<class T, class U> void SlaveReceiveGlobal(T& t, U& u)  {
        MPIBuffer buffer(MPISize(t) + MPISize(u));
        MPI_Bcast(buffer.GetBuffer(),buffer.GetSize(),MPI_DOUBLE,0,MPI_COMM_WORLD);
        buffer >> t >> u;
    }

    template<class T> void SlaveSendAdditive(const T& t) const  {
        MPIBuffer buffer(MPISize(t));
        buffer << t;
        MPI_Send(buffer.GetBuffer(),buffer.GetSize(),MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
    }

    template<class T> void MasterReceiveAdditive(T& t)  {

        for (int proc=1; proc<GetNprocs(); proc++)  {
            MPIBuffer buffer(MPISize(t));
            MPI_Status stat;
            MPI_Recv(buffer.GetBuffer(),buffer.GetSize(),MPI_DOUBLE,proc,TAG1,MPI_COMM_WORLD,&stat);
            t += buffer;
        }
    }

    template<class T> void MasterSendGeneArray(const ConstArray<T>& array) const    {

        int thusfar = 0;
        for (int proc=1; proc<GetNprocs(); proc++)  {
            int ngene = GetSlaveNgene(proc);
            MPIBuffer buffer(ngene* MPISize(array.GetVal(0)));
            for (int gene=0; gene<ngene; gene++)    {
                buffer << array.GetVal(thusfar);
                thusfar++;
            }
            MPI_Send(buffer.GetBuffer(),buffer.GetSize(),MPI_DOUBLE,proc,TAG1,MPI_COMM_WORLD);
        }
    }

    template<class T> void SlaveReceiveGeneArray(Array<T>& array)   {

        int ngene = GetLocalNgene();
        MPIBuffer buffer(ngene * MPISize(array[0]));
        MPI_Status stat;
        MPI_Recv(buffer.GetBuffer(),buffer.GetSize(),MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD,&stat);
        for (int gene=0; gene<ngene; gene++)    {
            buffer >> array[gene];
        }
    }

    template<class T> void SlaveSendGeneArray(const ConstArray<T>& array) const {

        int ngene = GetLocalNgene();
        MPIBuffer buffer(ngene * MPISize(array.GetVal(0)));
        for (int gene=0; gene<ngene; gene++)    {
            buffer << array.GetVal(gene);
        }
        MPI_Send(buffer.GetBuffer(),buffer.GetSize(),MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
    }

    template<class T> void MasterReceiveGeneArray(Array<T>& array)  {

        int thusfar = 0;
        for (int proc=1; proc<GetNprocs(); proc++)  {
            int ngene = GetSlaveNgene(proc);
            MPIBuffer buffer(ngene* MPISize(array[0]));
            MPI_Status stat;
            MPI_Recv(buffer.GetBuffer(),buffer.GetSize(),MPI_DOUBLE,proc,TAG1,MPI_COMM_WORLD,&stat);
            for (int gene=0; gene<ngene; gene++)    {
                buffer >> array[thusfar];
                thusfar++;
            }
        }
    }

    template<class T, class U> void MasterSendGeneArray(const ConstArray<T>& v, const ConstArray<U>& w) const    {

        int thusfar = 0;
        for (int proc=1; proc<GetNprocs(); proc++)  {
            int ngene = GetSlaveNgene(proc);
            MPIBuffer buffer(ngene* (MPISize(v.GetVal(0)) + MPISize(w.GetVal(0))));
            for (int gene=0; gene<ngene; gene++)    {
                buffer << v.GetVal(thusfar) << w.GetVal(thusfar);
                thusfar++;
            }
            MPI_Send(buffer.GetBuffer(),buffer.GetSize(),MPI_DOUBLE,proc,TAG1,MPI_COMM_WORLD);
        }
    }

    template<class T, class U> void SlaveReceiveGeneArray(Array<T>& v, Array<U>& w)   {

        int ngene = GetLocalNgene();
        MPIBuffer buffer(ngene * (MPISize(v[0]) + MPISize(w[0])));
        MPI_Status stat;
        MPI_Recv(buffer.GetBuffer(),buffer.GetSize(),MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD,&stat);
        for (int gene=0; gene<ngene; gene++)    {
            buffer >> v[gene] >> w[gene];
        }
    }

    template<class T, class U> void SlaveSendGeneArray(const ConstArray<T>& v, const ConstArray<U>& w) const {

        int ngene = GetLocalNgene();
        MPIBuffer buffer(ngene * (MPISize(v.GetVal(0)) + MPISize(w.GetVal(0))));
        for (int gene=0; gene<ngene; gene++)    {
            buffer << v.GetVal(gene) << w.GetVal(gene);
        }
        MPI_Send(buffer.GetBuffer(),buffer.GetSize(),MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
    }

    template<class T, class U> void MasterReceiveGeneArray(Array<T>& v, Array<U>& w)  {

        int thusfar = 0;
        for (int proc=1; proc<GetNprocs(); proc++)  {
            int ngene = GetSlaveNgene(proc);
            MPIBuffer buffer(ngene* (MPISize(v[0]) + MPISize(w[0])));
            MPI_Status stat;
            MPI_Recv(buffer.GetBuffer(),buffer.GetSize(),MPI_DOUBLE,proc,TAG1,MPI_COMM_WORLD,&stat);
            for (int gene=0; gene<ngene; gene++)    {
                buffer >> v[thusfar] >> w[thusfar];
                thusfar++;
            }
        }
    }

    protected:

    int myid;
    int nprocs;

	int Ngene;
    int LocalNgene;
    std::vector<int> SlaveNgene;
    std::vector<int> SlaveTotNsite;
    std::vector<int> GeneAlloc;
    std::vector<string> GeneName;
    std::vector<int> GeneNsite;

    SequenceAlignment* refdata;
};

#endif
