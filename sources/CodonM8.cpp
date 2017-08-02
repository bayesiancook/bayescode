#include <cmath>
#include <fstream>
#include "CodonM8Model.hpp"
using namespace std;

int main(int argc, char* argv[])	{

	string datafile = argv[1];
	string treefile = argv[2];
	int ncat = atoi(argv[3]);
	double pi = atof(argv[4]);
	string name = argv[5];

	CodonM8Model* model = new CodonM8Model(datafile,treefile,ncat,pi);
    model->Unfold();
	ofstream os((name + ".trace").c_str());
	ofstream pos((name + ".sitepp").c_str());
	model->TraceHeader(os);
	os.flush();
	while(1)	{
		model->Move();
		model->Trace(os);
        model->TracePostProb(pos);
		os.flush();
        pos.flush();
	}
}



