#ifndef CONTINUOUSDATA_H
#define CONTINUOUSDATA_H

#include "TaxonSet.hpp"
#include <fstream>
#include <cmath>
#include <sstream>

// this class works like an interface
// it does not do any job
class ContinuousData	{

	public:

	ContinuousData() {}

	ContinuousData(const TaxonSet* intaxset, int inNsite)	{
		taxset = intaxset;
		Nsite = inNsite;
		Data = new double*[GetNtaxa()];
		for (int i=0; i<GetNtaxa(); i++)	{
			Data[i] = new double[Nsite];
			for (int j=0; j<Nsite; j++)	{
				Data[i][j] = 0;
			}
		}
		charname.assign(GetNsite(), "none");
	}

	ContinuousData(ContinuousData* from)	{
		taxset = from->GetTaxonSet();
		Nsite = from->GetNsite();
		Data = new double*[GetNtaxa()];
		for (int i=0; i<GetNtaxa(); i++)	{
			Data[i] = new double[Nsite];
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
		// cerr << "??? in taxon set\n";
		os << GetNtaxa() << '\t' << GetNsite() << '\n';
		for (int i=0; i<GetNtaxa(); i++)	{
			os << GetTaxonSet()->GetTaxon(i);
			/*
			string s = taxset->GetTaxon(i);
			unsigned int l = s.length();
			unsigned int k = 0;
			while ((k < l) && (s[k] != '_')) k++;
			if (k == l)	{
				cerr << "error in get name\n";
				exit(1);
			}
			k++;
			os << s.substr(k,l-k);
			*/
			for (int j=0; j<GetNsite(); j++)	{
				os << '\t' << Data[i][j];
				/*
				if (Data[i][j] == -1)	{
					os << '\t' << -1;
				}
				else	{
					os << '\t' << log(Data[i][j]);
					// os << '\t' << exp(Data[i][j]);
				}
				*/
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

	void ToStreamLog(ostream& os) const {

		os << GetNtaxa() << '\t' << GetNsite() << '\n';
		for (int i=0; i<GetNtaxa(); i++)	{
			os << GetTaxonSet()->GetTaxon(i);
			/*
			string s = taxset->GetTaxon(i);
			unsigned int l = s.length();
			unsigned int k = 0;
			while ((k < l) && (s[k] != '_')) k++;
			if (k == l)	{
				cerr << "error in get name\n";
				exit(1);
			}
			k++;
			os << s.substr(k,l-k);
			*/
			for (int j=0; j<GetNsite(); j++)	{
				if (Data[i][j] == -1)	{
					os << '\t' << -1;
				}
				else	{
					os << '\t' << log(Data[i][j]);
				}
			}
			os << '\n';
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
	/*
		if (index == -1)	{
			cerr << "error : taxon not found : " << taxname << '\n';
			GetTaxonSet()->ToStream(cerr);
			exit(1);
		}
		return isMissing(index);
	}
	*/

	double GetMeanLog(int site) const {
		double total = 0;
		int n = 0;
		for (int j=0; j<GetNtaxa(); j++)	{
			if (Data[j][site] != -1)	{
				total += log(Data[j][site]);
				n++;
			}
		}
		return total / n;
	}

	string GetCharacterName(int site) const {
		return charname[site];
	}

	double GetState(string taxon, int site) const {
		return Data[taxset->GetTaxonIndex(taxon)][site];
	}

	double GetState(int taxon, int site) const {
		return Data[taxon][site];
	}

	int Nsite;
	const TaxonSet* taxset;
	double** Data;
	vector<string> charname;

};

class FileContinuousData : public ContinuousData {

	public:
		FileContinuousData(istream& is)	{
			ReadDataFromFile(is);
		}

		FileContinuousData(string filename)	{
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
		Data = new double*[Ntaxa];
		for (int i=0; i<Ntaxa; i++)	{
			is >> name[i];
			Data[i] = new double[Nsite];
			for (int j=0; j<Nsite; j++)	{
				is >> Data[i][j];
			}
		}
		taxset = new TaxonSet(name);
		return 1;
	}
};

#endif
