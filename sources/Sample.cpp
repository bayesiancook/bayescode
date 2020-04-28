
#include "Sample.hpp"

Sample::Sample(string filename, int in_burnin, int in_every, int in_until) {
    burnin = in_burnin;
    every = in_every;
    until = in_until;
    name = filename;
    chain_is = 0;
    chainsaveall = 1;
}

Sample::~Sample() { delete chain_is; }

void Sample::OpenChainFile() {
    if (!chainsaveall) {
        cerr << "error: no chain file\n";
        exit(1);
    }
    if (until == -1) {
        until = chainsize;
    }
    if (until > chainsize) {
        cerr << "number of points saved is less than " << until << '\n';
        until = chainsize;
    }
    if (burnin == -1) {
        burnin = chainsize / 10;
    }
    size = (until - burnin) / every;
    if (size <= 0) {
        cerr << "error : chain not long enough\n";
        cerr << burnin << '\t' << every << '\t' << until << '\n';
        exit(1);
    }
    currentpoint = 0;

    chain_is = new ifstream((name + ".chain").c_str());
    if (!chain_is) {
        cerr << "error: cannot find file " << name << ".chain\n";
        exit(1);
    }
    string line;
    getline(*chain_is,line);
    for (int i = 0; i < burnin; i++) {
        model->FromStream(*chain_is);
    }
}

void Sample::GetNextPoint() {
    if (currentpoint == size) {
        cerr << "error in Sample::GetNextPoint: going past last points\n";
        exit(1);
    }
    if (currentpoint) {
        for (int i = 0; i < every - 1; i++) {
            model->FromStream(*chain_is);
        }
    }
    model->FromStream(*chain_is);
    currentpoint++;
}

void Sample::PostPred() {
    cerr << size << " points to read\n";
    for (int i = 0; i < size; i++) {
        cerr << '.';
        GetNextPoint();
        ostringstream s;
        s << "ppred" << name << "_" << i;
        // s << "ppred" << name << "_" << i << ".ali";
        model->PostPred(s.str());
    }
    cerr << '\n';
}

void Sample::AllPostPred() {

    vector<string> statnames = GetModel()->GetPostPredStatNames();
    size_t nstat = statnames.size();
    vector<double> stats(nstat, 0);
    vector<double> obsstats(nstat, 0);
    vector<double> meanstats(nstat, 0);
    vector<double> varstats(nstat, 0);
    vector<double> ppstats(nstat, 0);

    model->AllPost(obsstats);

    cerr << size << " points to read\n";
    for (int i=0; i<size; i++) {
        cerr << '.';
        GetNextPoint();

        model->AllPostPred(stats);
        for (size_t j=0; j<nstat; j++)  {
            meanstats[j] += stats[j];
            varstats[j] += stats[j]*stats[j];
            if (stats[j] > obsstats[j]) {
                ppstats[j]++;
            }
        }
    }
    cerr << '\n';

    ofstream os((name + ".ppred").c_str());
    os << "stat\tobs\tpred\tzscore\tpp\n";
    for (size_t j=0; j<nstat; j++)  {
        ppstats[j] /= size;
        meanstats[j] /= size;
        varstats[j] /= size;
        varstats[j] -= meanstats[j]*meanstats[j];
        double z = (obsstats[j] - meanstats[j]) / sqrt(varstats[j]);
        os << statnames[j] << '\t' << obsstats[j] << '\t' << meanstats[j] << '\t' << z << '\t' << ppstats[j] << '\n';
    }
}

