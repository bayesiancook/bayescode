
#include "CodonSequenceAlignment.hpp"
#include "Random.hpp"

int main(int argc, char* argv[])	{

    if (argc == 1)  {
        cerr << "ppredtest obs_ali basename nrep\n";
        exit(1);
    }
	string obsali = argv[1];
	string basename = argv[2];
	int nrep = atoi(argv[3]);

	FileSequenceAlignment ali(obsali);
	CodonSequenceAlignment codali(&ali);
	double obsdiv = codali.GetMeanAADiversity();
	double obsdiff = codali.GetMeanDiff();
	double obsdnds = codali.GetMeanEmpiricaldNdS();

	double meandiv = 0;
	double vardiv = 0;
	double ppdiv = 0;

	double meandiff = 0;
	double vardiff = 0;
	double ppdiff = 0;

	double meandnds = 0;
	double vardnds = 0;
	double ppdnds = 0;

	for (int rep=0; rep<nrep; rep++)	{
		cerr << '.';
		ostringstream s;
		s << "ppred" << basename << "_" << rep << ".ali";
		FileSequenceAlignment ali(s.str());
		CodonSequenceAlignment codali(&ali);
		double div = codali.GetMeanAADiversity();
		meandiv += div;
		vardiv += div*div;
		if (div < obsdiv)	{
			ppdiv++;
		}

		double diff = codali.GetMeanDiff();
		meandiff += diff;
		vardiff += diff*diff;
		if (diff > obsdiff)	{
			ppdiff++;
		}

		double dnds = codali.GetMeanEmpiricaldNdS();
		meandnds += dnds;
		vardnds += dnds*dnds;
		if (dnds > obsdnds)	{
			ppdnds++;
		}
	}
	cerr << '\n';

	meandiv /= nrep;
	vardiv /= nrep;
	vardiv -= meandiv*meandiv;
	ppdiv /= nrep;
	double zdiv = (meandiv - obsdiv) / sqrt(vardiv);
	cout << '\n';
	cout << "mean site-specific amino-acid diversity\n";
	cout << "obs  : " << obsdiv << '\n';
	cout << "pred : " << meandiv << '\n';
	cout << "z    : " << zdiv << '\n';
	cout << "pp   : " << ppdiv << '\n';

	cout << '\n';
	cout << "mean pairwise divergence\n";
	meandiff /= nrep;
	vardiff /= nrep;
	vardiff -= meandiff*meandiff;
	ppdiff /= nrep;
	double zdiff = (meandiff - obsdiff) / sqrt(vardiff);
	cout << "obs  : " << obsdiff << '\n';
	cout << "pred : " << meandiff << '\n';
	cout << "z    : " << zdiff << '\n';
	cout << "pp   : " << ppdiff << '\n';

	cout << '\n';
	cout << "mean empirical dNdS\n";
	meandnds /= nrep;
	vardnds /= nrep;
	vardnds -= meandnds*meandnds;
	ppdnds /= nrep;
	double zdnds = (meandnds - obsdnds) / sqrt(vardnds);
	cout << "obs  : " << obsdnds << '\n';
	cout << "pred : " << meandnds << '\n';
	cout << "z    : " << zdnds << '\n';
	cout << "pp   : " << ppdnds << '\n';
	cout << '\n';

}


