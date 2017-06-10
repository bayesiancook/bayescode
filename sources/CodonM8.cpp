#include <cmath>
#include <fstream>
#include "CodonM8Model.hpp"
using namespace std;

int main(int argc, char* argv[])	{

	string datafile = argv[1];
	string treefile = argv[2];
    int withpos = atoi(argv[3]);
	string name = argv[4];

	CodonM8Model* model = new CodonM8Model(datafile,treefile,withpos);
	ofstream os((name + ".trace").c_str());
	model->TraceHeader(os);
	os.flush();
	while(1)	{
		model->Move();
		model->Trace(os);
		os.flush();
	}
}



