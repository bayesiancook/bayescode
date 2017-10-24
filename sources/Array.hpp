
#ifndef ARRAY_H
#define ARRAY_H

#include <vector>
#include "MPIBuffer.hpp"

template<class T> class ConstArray	{

	public:
	virtual ~ConstArray() {}

	virtual int GetSize() const = 0;
	virtual const T& GetVal(int index) const = 0;

    unsigned int GetMPISize() const {return this->GetSize() * MPISize(this->GetVal(0));}

    void MPIPut(MPIBuffer& buffer) const {
        for (int i=0; i<this->GetSize(); i++)  {
            buffer << this->GetVal(i);
        }
    }

    void ToStream(ostream& os) const    {
        for (int i=0; i<this->GetSize(); i++)  {
            os << this->GetVal(i) << '\t';;
        }
	os << '\n';
    }
};

template<class T> class Array : public ConstArray<T>	{

	public:
	virtual ~Array() {}

	virtual T& operator[](int index) = 0;

    void MPIGet(const MPIBuffer& buffer)    {
        for (int i=0; i<this->GetSize(); i++)  {
            buffer >> (*this)[i];
        }
    }

    void FromStream(istream& is)    {
        for (int i=0; i<this->GetSize(); i++)  {
            is >> (*this)[i];
        }
    }
};

template<class T> ostream& operator<<(ostream& os, const ConstArray<T>& array)  {
    array.ToStream(os);
    return os;
}

template<class T> istream& operator>>(istream& is, Array<T>& array) {
    array.FromStream(is);
    return is;
}

template<class T> ostream& operator<<(ostream& os, const vector<T>& array)  {
    for (unsigned int i=0; i<array.size(); i++) {
        os << array[i] << '\t';
    }
    os << '\n';
    return os;
}

template<class T> istream& operator>>(istream& is, vector<T>& array)    {
    for (unsigned int i=0; i<array.size(); i++) {
        is >> array[i];
    }
    return is;
}

template<class T> class HomogeneousArray : public ConstArray<T>	{

	public:
	HomogeneousArray(int insize, const T& invalue) : size(insize), value(invalue) {}
	~HomogeneousArray() {}

	int GetSize() const /*override*/ {return size;}
	const T& GetVal(int index) const /*override*/ {return value;}

	private:
	int size;
	const T& value;
};

template<class T> class SimpleArray : public Array<T>	{

	public:
	SimpleArray(int insize) : array(insize) {}
    SimpleArray(int insize, const T& initval) : array(insize,initval) {}
	virtual ~SimpleArray() {}

	int GetSize() const /*override*/ {return array.size();}
	T& operator[](int index) /*override*/ {return array[index];}
	const T& GetVal(int index) const /*override*/ {return array[index];}
	const vector<T>& GetArray() const {return array;}

    void Swap(int cat1, int cat2)   {
        T tmp = (*this)[cat1];
        (*this)[cat1] = (*this)[cat2];
        (*this)[cat2] = tmp;
    }

	protected:
	vector<T> array;
};


template<class T> class ConstMixtureArray : public ConstArray<T>	{

	public:

	ConstMixtureArray(const Array<T>* incomponents, const Array<int>* inalloc) : components(incomponents), alloc(inalloc)	{
	}
	~ConstMixtureArray() {}

	int GetSize() const /*override*/ {return alloc->GetSize();}
	const T& GetVal(int i) const /*override*/ {
		return components->GetVal(alloc->GetVal(i));
	}

	private:
	const Array<T>* components;
	const Array<int>* alloc;
};

template<class T> class MixtureArray : public Array<T>	{

	public:

	MixtureArray(Array<T>* incomponents, const Array<int>* inalloc) : components(incomponents), alloc(inalloc)	{
	}
	~MixtureArray() {}

	int GetSize() const /*override*/ {return alloc->GetSize();}
	const T& GetVal(int i) const /*override*/ {
		return components->GetVal(alloc->GetVal(i));
	}
	T& operator[](int i) /*override*/ {return (*components)[alloc->GetVal(i)];}

	private:
	const Array<T>* components;
	const Array<int>* alloc;
};

#endif

