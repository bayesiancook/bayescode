#ifndef INTEGERDATA_H
#define INTEGERDATA_H

#include "TaxonSet.hpp"
#include <fstream>
#include <cmath>
#include <sstream>

// this class works like an interface
// it does not do any job
class IntegerData	{

	public:

	IntegerData() {}

	IntegerData(const TaxonSet* intaxset, int inNsite)	{
		taxset = intaxset;
		Nsite = inNsite;
		Data = new int*[GetNtaxa()];
		for (int i=0; i<GetNtaxa(); i++)	{
			Data[i] = new int[Nsite];
			for (int j=0; j<Nsite; j++)	{
				Data[i][j] = -1;
			}
		}
		charname.assign(GetNsite(), "none");
	}

	IntegerData(IntegerData* from)	{
		taxset = from->GetTaxonSet();
		Nsite = from->GetNsite();
		Data = new int*[GetNtaxa()];
		for (int i=0; i<GetNtaxa(); i++)	{
			Data[i] = new int[Nsite];
			for (int j=0; j<Nsite; j++)	{
				Data[i][j] = from->Data[i][j];
			}
		}
		charname.assign(GetNsite(), "");
		for (int j=0; j<Nsite; j++)	{
			charname[j] = from->charname[j];
		}
	}

	void ToStream(ostream& os, TaxonSet* taxset = 0) const {

		if (! taxset)	{
            os << GetNtaxa() << '\t' << GetNsite() << '\n';
            for (int i=0; i<GetNtaxa(); i++)	{
                os << GetTaxonSet()->GetTaxon(i);
                for (int j=0; j<GetNsite(); j++)	{
                    os << '\t' << Data[i][j];
                }
                os << '\n';
            }
		}
		else	{
			int ntaxa = 0;
			for (int i=0; i<GetNtaxa(); i++)	{
				if (taxset->GetTaxonIndex(GetTaxonSet()->GetTaxon(i)) != -1)	{
					ntaxa++;
				}
			}

			os << ntaxa << '\t' << GetNsite() << '\n';
			for (int i=0; i<GetNtaxa(); i++)	{
				if (taxset->GetTaxonIndex(GetTaxonSet()->GetTaxon(i)) != -1)	{
					os << GetTaxonSet()->GetTaxon(i);
					for (int j=0; j<GetNsite(); j++)	{
						os << '\t' << Data[i][j];
					}
					os << '\n';
				}
			}
		}
	}

	// the list of taxa
	const TaxonSet* GetTaxonSet() const {
		return taxset;
	}

	int GetNtaxa() const {
		return GetTaxonSet()->GetNtaxa();
	}

	int GetNsite() const {
		return Nsite;
	}

	bool isMissing(int taxon, int site) const {
		return Data[taxon][site] == -1;
	}

	bool isMissing(int taxon) const {
		bool mis = true;
		for (int i=0; i<Nsite; i++)	{
			mis &= (Data[taxon][i] == -1);
		}
		return mis;
	}

	bool isMissing(string taxname) const {
		int index = GetTaxonSet()->GetTaxonIndex(taxname);
		return (index == -1);
	}

	string GetCharacterName(int site) const {
		return charname[site];
	}

	int GetState(string taxon, int site) const {
		return Data[taxset->GetTaxonIndex(taxon)][site];
	}

	int GetState(int taxon, int site) const {
		return Data[taxon][site];
	}

	int Nsite;
	const TaxonSet* taxset;
	int** Data;
	vector<string> charname;

};

class FileIntegerData : public IntegerData {

	public:
		FileIntegerData(istream& is)	{
			ReadDataFromFile(is);
		}

		FileIntegerData(string filename)	{
			ifstream is(filename.c_str());
			if (! is)	{
				cerr << "error when opening file : " << filename << '\n';
				exit(1);
			}
			ReadDataFromFile(is);
		}

	private:

	int  ReadDataFromFile(istream& is)	{

		string temp;
		is >> temp;
		int Ntaxa;
		if (temp == "#TRAITS")	{
			is >> Ntaxa;
			is >> Nsite;
			charname.assign(Nsite, "");
			for (int j=0; j<Nsite; j++)	{
				is >> charname[j];
			}
		}
		else	{
			Ntaxa = atoi(temp.c_str());
			is >> Nsite;
			charname.assign(Nsite, "");
			for (int j=0; j<Nsite; j++)	{
				ostringstream s;
				s << "character" << j+1;
				charname[j] = s.str();
			}
		}


		vector<string> name(Ntaxa, "");
		Data = new int*[Ntaxa];
		for (int i=0; i<Ntaxa; i++)	{
			is >> name[i];
			Data[i] = new int[Nsite];
			for (int j=0; j<Nsite; j++)	{
				double tmp;
				is >> tmp;
				Data[i][j] = int(tmp);
			}
		}
		taxset = new TaxonSet(name);
		return 1;
	}
};

#endif
