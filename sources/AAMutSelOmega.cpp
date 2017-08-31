#include <cmath>
#include <fstream>
#include "AAMutSelOmegaModel.hpp"
using namespace std;

int main(int argc, char* argv[])	{

	string datafile = argv[1];
	string treefile = argv[2];
	string name = argv[3];

	AAMutSelOmegaModel* model = new AAMutSelOmegaModel(datafile,treefile);
    model->Allocate();
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


