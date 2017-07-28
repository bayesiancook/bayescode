#include <cmath>
#include <fstream>
#include "SingleOmegaModel.hpp"
using namespace std;

int main(int argc, char* argv[])	{

	string datafile = argv[1];
	string treefile = argv[2];
	int withnuc = atoi(argv[3]);
	string name = argv[4];

	SingleOmegaModel* model = new SingleOmegaModel(datafile,treefile,withnuc);
	model->Unfold();
	ofstream os((name + ".trace").c_str());
	model->TraceHeader(os);
	
	os.flush();
	while(1)	{
		model->Move();
		model->Trace(os);
		os.flush();
	}
}


