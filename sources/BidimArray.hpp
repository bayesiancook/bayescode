
#ifndef BIDIMARRAY_H
#define BIDIMARRAY_H

#include <vector>
#include "MPIBuffer.hpp"

template<class T> class ConstBidimArray	{

	public:
	virtual ~ConstBidimArray() {}

	virtual int GetNrow() const = 0;
    virtual int GetNcol() const = 0;
	virtual const T& GetVal(int i, int j) const = 0;

    unsigned int GetMPISize() const {return this->GetNrow() * this->GetNcol() * MPISize(this->GetVal(0));}

    void MPIPut(MPIBuffer& buffer) const {
        for (int i=0; i<this->GetNrow(); i++)  {
            for (int j=0; j<this->GetNcol(); j++)   {
                buffer << this->GetVal(i,j);
            }
        }
    }

    void ToStream(ostream& os) const    {
        for (int i=0; i<this->GetNrow(); i++)  {
            for (int j=0; j<this->GetNcol(); j++)   {
                os << this->GetVal(i,j);
            }
        }
    }
};

template<class T> class BidimArray : public ConstBidimArray<T>	{

	public:
	virtual ~BidimArray() {}

	virtual T& operator()(int i, int j) = 0;

    void MPIGet(const MPIBuffer& buffer)    {
        for (int i=0; i<this->GetNrow(); i++)  {
            for (int j=0; j<this->GetNcol(); j++)   {
                buffer >> (*this)(i,j);
            }
        }
    }

    void FromStream(istream& is)    {
        for (int i=0; i<this->GetNrow(); i++)  {
            for (int j=0; j<this->GetNcol(); j++)   {
                is >> (*this)(i,j);
            }
        }
    }
};

template<class T> ostream& operator<<(ostream& os, const ConstBidimArray<T>& array)  {
    array.ToStream(os);
    return os;
}

template<class T> istream& operator>>(istream& is, BidimArray<T>& array) {
    array.FromStream(is);
    return is;
}

template<class T> class HomogeneousBidimArray : public ConstBidimArray<T>	{

	public:
	HomogeneousBidimArray(int innrow, int inncol, const T& invalue) : nrow(innrow), ncol(inncol), value(invalue) {}
	~HomogeneousBidimArray() {}

    int GetNrow() const override {return nrow;}
    int GetNcol() const override {return ncol;}

	const T& GetVal(int i, int j) const override {return value;}

	private:
    int nrow;
    int ncol;
	const T& value;
};

template<class T> class SimpleBidimArray : public BidimArray<T>	{

	public:
	SimpleBidimArray(int innrow, int inncol, const T& initval) : nrow(innrow), ncol(inncol), array(innrow,vector<T>(inncol,initval)) {}
	virtual ~SimpleBidimArray() {}

    int GetNrow() const override {return nrow;}
    int GetNcol() const override {return ncol;}

	T& operator()(int i, int j) override {return array[i][j];}
	const T& GetVal(int i, int j) const override {return array[i][j];}

	const vector<T>& GetSubArray(int i) const {return array[i];}

	protected:
    int nrow;
    int ncol;
	vector<vector<T> > array;
};

#endif

