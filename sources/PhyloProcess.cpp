#include "PhyloProcess.hpp"
#include <algorithm>
#include "PathSuffStat.hpp"
using namespace std;

PhyloProcess::PhyloProcess(const Tree *intree, const SequenceAlignment *indata,
                           const BranchSelector<double> *inbranchlength,
                           const Selector<double> *insiterate,
                           const BranchSiteSelector<SubMatrix> *insubmatrixarray,
                           const Selector<SubMatrix> *inrootsubmatrixarray) {
    tree = intree;
    data = indata;
    Nstate = data->GetNstate();
    maxtrial = DEFAULTMAXTRIAL;
    branchlength = inbranchlength;
    siterate = insiterate;
    submatrixarray = insubmatrixarray;
    allocsubmatrixarray = false;
    rootsubmatrixarray = inrootsubmatrixarray;
    allocrootsubmatrixarray = false;
}

PhyloProcess::PhyloProcess(const Tree *intree, const SequenceAlignment *indata,
                           const BranchSelector<double> *inbranchlength,
                           const Selector<double> *insiterate, const SubMatrix *insubmatrix) {
    tree = intree;
    data = indata;
    Nstate = data->GetNstate();
    maxtrial = DEFAULTMAXTRIAL;
    branchlength = inbranchlength;
    siterate = insiterate;
    submatrixarray =
        new BranchHomogeneousSiteHomogeneousSelector<SubMatrix>(*tree, GetNsite(), *insubmatrix);
    allocsubmatrixarray = true;
    rootsubmatrixarray = new HomogeneousSelector<SubMatrix>(GetNsite(), *insubmatrix);
    allocrootsubmatrixarray = true;
}

PhyloProcess::PhyloProcess(const Tree *intree, const SequenceAlignment *indata,
                           const BranchSelector<double> *inbranchlength,
                           const Selector<double> *insiterate,
                           const Selector<SubMatrix> *insubmatrixarray) {
    tree = intree;
    data = indata;
    Nstate = data->GetNstate();
    maxtrial = DEFAULTMAXTRIAL;
    branchlength = inbranchlength;
    siterate = insiterate;
    if (insubmatrixarray->GetSize() != GetNsite()) {
        std::cerr << "error in PhyloProcess constructor: size of matrix array does "
                     "not match alignment size\n";
        exit(1);
    }
    submatrixarray =
        new BranchHomogeneousSiteHeterogeneousSelector<SubMatrix>(*tree, *insubmatrixarray);
    allocsubmatrixarray = true;
    rootsubmatrixarray = insubmatrixarray;
    allocrootsubmatrixarray = false;
}

PhyloProcess::PhyloProcess(const Tree *intree, const SequenceAlignment *indata,
                           const BranchSelector<double> *inbranchlength,
                           const Selector<double> *insiterate,
                           const BranchSelector<SubMatrix> *insubmatrixbrancharray,
                           const SubMatrix *insubmatrix) {
    tree = intree;
    data = indata;
    Nstate = data->GetNstate();
    maxtrial = DEFAULTMAXTRIAL;
    branchlength = inbranchlength;
    siterate = insiterate;
    submatrixarray = new BranchHeterogeneousSiteHomogeneousSelector<SubMatrix>(
        *insubmatrixbrancharray, GetNsite());
    allocsubmatrixarray = true;
    rootsubmatrixarray = new HomogeneousSelector<SubMatrix>(GetNsite(), *insubmatrix);
    allocrootsubmatrixarray = true;
}

PhyloProcess::~PhyloProcess() {
    Cleanup();
    if (allocsubmatrixarray) {
        delete submatrixarray;
    }
    if (allocrootsubmatrixarray) {
        delete rootsubmatrixarray;
    }
}

void PhyloProcess::SetData(const SequenceAlignment *indata) { data = indata; }

void PhyloProcess::Unfold() {
    sitearray = new int[GetNsite()];
    sitelnL = new double[GetNsite()];
    for (int i = 0; i < GetNsite(); i++) {
        sitearray[i] = 1;
    }
    CreateMissingMap();
    FillMissingMap();
    RecursiveCreate(GetRoot());
    RecursiveCreateTBL(GetRoot());
    ClampData();
}

void PhyloProcess::Cleanup() {
    DeleteMissingMap();
    RecursiveDeleteTBL(GetRoot());
    RecursiveDelete(GetRoot());
    delete[] sitearray;
    delete[] sitelnL;
}

void PhyloProcess::CreateMissingMap() {
    missingmap = new int *[GetTree()->GetNnode()];
    for (int j = 0; j < GetTree()->GetNnode(); j++) {
        missingmap[j] = new int[GetNsite()];
        for (int i = 0; i < GetNsite(); i++) {
            missingmap[j][i] = -1;
        }
    }
}

void PhyloProcess::DeleteMissingMap() {
    for (int j = 0; j < GetTree()->GetNnode(); j++) {
        delete[] missingmap[j];
    }
    delete[] missingmap;
}

void PhyloProcess::FillMissingMap() {
    /*
        for (int j=0; j<GetTree()->GetNnode(); j++)	{
                for (int i=0; i<GetNsite(); i++)	{
                        missingmap[j][i] = 0;
                }
        }
    int rootindex = GetTree()->GetRoot()->GetNode()->GetIndex();
    for (int i=0; i<GetNsite(); i++)	{
        missingmap[rootindex][i] = 2;
    }
    */

    BackwardFillMissingMap(GetRoot());
    ForwardFillMissingMap(GetRoot(), GetRoot());
}

/*
void PhyloProcess::PrintMissing(const Link* from, int site)    {

    if (!from->isLeaf()) {
        cerr << "(";
        for (const Link* link=from->Next(); link!=from; link=link->Next())  {
            PrintMissing(link->Out(),site);
            if (link->Next() != from)   {
                cerr << ",";
            }
        }
        cerr << ")";
    }
    cerr << missingmap[from->GetNode()->GetIndex()][site];
}
*/

void PhyloProcess::BackwardFillMissingMap(const Link *from) {
    int index = from->GetNode()->GetIndex();
    for (int i = 0; i < GetNsite(); i++) {
        missingmap[index][i] = 0;
    }
    if (from->isLeaf()) {
        for (int i = 0; i < GetNsite(); i++) {
            int state = GetData(index, i);
            if (state != -1) {
                missingmap[index][i] = 1;
            }
        }
    } else {
        for (const Link *link = from->Next(); link != from; link = link->Next()) {
            BackwardFillMissingMap(link->Out());
            int j = link->Out()->GetNode()->GetIndex();
            for (int i = 0; i < GetNsite(); i++) {
                if (missingmap[j][i]) {
                    missingmap[index][i]++;
                }
            }
        }
    }
}

void PhyloProcess::ForwardFillMissingMap(const Link *from, const Link *up) {
    int index = from->GetNode()->GetIndex();
    int upindex = up->GetNode()->GetIndex();

    if (from->isRoot()) {
        for (int i = 0; i < GetNsite(); i++) {
            if (missingmap[index][i] <= 1) {
                missingmap[index][i] = 0;
            } else {
                missingmap[index][i] = 2;
            }
        }
    } else {
        for (int i = 0; i < GetNsite(); i++) {
            if (missingmap[index][i] > 0) {
                if (missingmap[upindex][i]) {
                    missingmap[index][i] = 1;
                } else {
                    if (from->isLeaf() || (missingmap[index][i] > 1)) {
                        missingmap[index][i] = 2;
                    } else {
                        missingmap[index][i] = 0;
                    }
                }
            }
        }
    }
    for (const Link *link = from->Next(); link != from; link = link->Next()) {
        ForwardFillMissingMap(link->Out(), from);
    }
}

void PhyloProcess::RecursiveCreate(const Link *from) {
    auto state = new int[GetNsite()];
    statemap[from->GetNode()] = state;

    auto array = new BranchSitePath *[GetNsite()];
    for (int i = 0; i < GetNsite(); i++) {
        array[i] = 0;
    }
    pathmap[from->GetNode()] = array;

    for (const Link *link = from->Next(); link != from; link = link->Next()) {
        RecursiveCreate(link->Out());
    }
}

void PhyloProcess::RecursiveDelete(const Link *from) {
    for (const Link *link = from->Next(); link != from; link = link->Next()) {
        RecursiveDelete(link->Out());
    }

    delete[] statemap[from->GetNode()];

    BranchSitePath **path = pathmap[from->GetNode()];
    for (int i = 0; i < GetNsite(); i++) {
        delete path[i];
    }
}

void PhyloProcess::RecursiveCreateTBL(const Link *from) {
    condlmap[from] = new double[GetNstate() + 1];
    for (Link *link = from->Next(); link != from; link = link->Next()) {
        condlmap[link] = new double[GetNstate() + 1];
        RecursiveCreateTBL(link->Out());
    }
}

void PhyloProcess::RecursiveDeleteTBL(const Link *from) {
    for (Link *link = from->Next(); link != from; link = link->Next()) {
        RecursiveDeleteTBL(link->Out());
        delete[] condlmap[link];
    }
    delete[] condlmap[from];
}

double PhyloProcess::SiteLogLikelihood(int site) const {
    Pruning(GetRoot(), site);
    double ret = 0;
    double *t = GetCondLikelihood(GetRoot());
    const EVector &stat = GetRootFreq(site);

    for (int k = 0; k < GetNstate(); k++) {
        ret += t[k] * stat[k];
    }
    if (ret == 0) {
        cerr << "pruning : 0 \n";
        for (int k = 0; k < GetNstate(); k++) {
            cerr << t[k] << '\t' << stat[k] << '\n';
        }
        exit(1);
    }
    return log(ret) + t[GetNstate()];
}

double PhyloProcess::FastSiteLogLikelihood(int site) const {
    double ret = 0;
    double *t = GetCondLikelihood(GetRoot());
    const EVector &stat = GetRootFreq(site);
    for (int k = 0; k < GetNstate(); k++) {
        ret += t[k] * stat[k];
    }
    if (ret == 0) {
        cerr << "pruning : 0 \n";
        for (int k = 0; k < GetNstate(); k++) {
            cerr << t[k] << '\t' << stat[k] << '\n';
        }
        exit(1);
    }
    sitelnL[site] = log(ret) + t[GetNstate()];
    return sitelnL[site];
}

double PhyloProcess::GetFastLogProb() const {
    double total = 0;
    MeasureTime timer;
    for (int i = 0; i < GetNsite(); i++) {
        total += sitelnL[i];
    }
    // timer.print<2>("GetFastLogProb. ");
    return total;
}

double PhyloProcess::GetLogLikelihood() const {
#if DEBUG > 1
    MeasureTime timer;
#endif
    double total = 0;
    for (int i = 0; i < GetNsite(); i++) {
        total += SiteLogLikelihood(i);
    }
#if DEBUG > 1
    timer.print<2>("GetLogProb. ");
#endif
    return total;
}

void PhyloProcess::Pruning(const Link *from, int site) const {
    double *t = GetCondLikelihood(from);
    if (from->isLeaf()) {
        int totcomp = 0;
        for (int k = 0; k < GetNstate(); k++) {
            // if (polyprocess != nullptr) {
            //    t[k] = polyprocess->GetLogProb(site,from->GetNode()->GetIndex(),k);
            //}
            // else    {
            if (isDataCompatible(from->GetNode()->GetIndex(), site, k)) {
                t[k] = 1;
                totcomp++;
            } else {
                t[k] = 0;
            }
            //}
        }
        if (totcomp == 0) {
            cerr << "error : no compatibility\n";
            cerr << GetData(from->GetNode()->GetIndex(), site) << '\n';
            exit(1);
        }

        t[GetNstate()] = 0;
    } else {
        for (int k = 0; k < GetNstate(); k++) {
            t[k] = 1.0;
        }
        t[GetNstate()] = 0;
        for (const Link *link = from->Next(); link != from; link = link->Next()) {
            double *tbl = GetCondLikelihood(link);
            Pruning(link->Out(), site);
            GetSubMatrix(link->GetBranch()->GetIndex(), site)
                .BackwardPropagate(
                    GetCondLikelihood(link->Out()), GetCondLikelihood(link),
                    GetBranchLength(link->GetBranch()->GetIndex()) * GetSiteRate(site));
            for (int k = 0; k < GetNstate(); k++) {
                t[k] *= tbl[k];
            }
            t[GetNstate()] += tbl[GetNstate()];
        }
        double max = 0;
        for (int k = 0; k < GetNstate(); k++) {
            if (t[k] < 0) {
                /*
                  cerr << "error in pruning: negative prob : " << t[k] << "\n";
                  exit(1);
                */
                t[k] = 0;
            }
            if (max < t[k]) {
                max = t[k];
            }
        }
        if (max == 0) {
            cerr << "max = 0\n";
            cerr << "error in pruning: null likelihood\n";
            if (from->isRoot()) {
                cerr << "is root\n";
            }
            cerr << '\n';
            exit(1);
            max = 1e-20;
        }
        for (int k = 0; k < GetNstate(); k++) {
            t[k] /= max;
        }
        t[GetNstate()] += log(max);
    }
}

void PhyloProcess::PruningAncestral(const Link *from, int site) {
    int &state = GetState(from->GetNode(), site);
    if (from->isRoot()) {
        auto aux = new double[GetNstate()];
        auto cumulaux = new double[GetNstate()];
        try {
            double *tbl = GetCondLikelihood(from);
            const EVector &stat = GetRootFreq(site);
            double tot = 0;
            for (int k = 0; k < GetNstate(); k++) {
                aux[k] = stat[k] * tbl[k];
                tot += aux[k];
                cumulaux[k] = tot;
            }
            double u = tot * Random::Uniform();
            int s = 0;
            while ((s < GetNstate()) && (cumulaux[s] < u)) {
                s++;
            }
            if (s == GetNstate()) {
                cerr << "error in pruning ancestral: overflow\n";
                exit(1);
            }
            state = s;
        } catch (...) {
            cerr << "in root::PruningAncestral\n";
            for (int k = 0; k < GetNstate(); k++) {
                cerr << aux[k] << '\n';
            }
            exit(1);
            throw;
        }
        delete[] aux;
        delete[] cumulaux;
    }
    for (const Link *link = from->Next(); link != from; link = link->Next()) {
        auto aux = new double[GetNstate()];
        auto cumulaux = new double[GetNstate()];
        try {
            for (int k = 0; k < GetNstate(); k++) {
                aux[k] = 1;
            }
            GetSubMatrix(link->GetBranch()->GetIndex(), site)
                .GetFiniteTimeTransitionProb(
                    state, aux, GetBranchLength(link->GetBranch()->GetIndex()) * GetSiteRate(site));
            double *tbl = GetCondLikelihood(link->Out());
            for (int k = 0; k < GetNstate(); k++) {
                aux[k] *= tbl[k];
            }

            // dealing with numerical problems:
            double max = 0;
            for (int k = 0; k < GetNstate(); k++) {
                if (aux[k] < 0) {
                    aux[k] = 0;
                }
                if (max < aux[k]) {
                    max = aux[k];
                }
            }
            if (max == 0) {
                auto stat = GetSubMatrix(link->GetBranch()->GetIndex(), site).GetStationary();
                for (int k = 0; k < GetNstate(); k++) {
                    aux[k] = stat[k];
                }
            }
            // end of dealing with dirty numerical problems
            double tot = 0;
            for (int k = 0; k < GetNstate(); k++) {
                tot += aux[k];
                cumulaux[k] = tot;
            }
            double u = tot * Random::Uniform();
            int s = 0;
            while ((s < GetNstate()) && (cumulaux[s] < u)) {
                s++;
            }
            if (s == GetNstate()) {
                cerr << "error in pruning ancestral: overflow\n";
                exit(1);
            }
            int &nodestate = GetState(link->Out()->GetNode(), site);
            nodestate = s;
        } catch (...) {
            cerr << "in internal leave::PruningAncestral\n";
            for (int k = 0; k < GetNstate(); k++) {
                cerr << aux[k] << '\n';
            }
            exit(1);
            throw;
        }
        delete[] aux;
        delete[] cumulaux;
        PruningAncestral(link->Out(), site);
    }
}

void PhyloProcess::RootPosteriorDraw(int site) {
    auto aux = new double[GetNstate()];
    double *tbl = GetCondLikelihood(GetRoot());
    const EVector &stat = GetRootFreq(site);
    for (int k = 0; k < GetNstate(); k++) {
        aux[k] = stat[k] * tbl[k];
    }
    GetState(GetRoot()->GetNode(), site) = Random::DrawFromDiscreteDistribution(aux, GetNstate());
    delete[] aux;
}

void PhyloProcess::PriorSample(const Link *from, int site, bool rootprior) {
    int &state = GetState(from->GetNode(), site);
    if (from->isRoot()) {
        if (rootprior) {
            state = Random::DrawFromDiscreteDistribution(GetRootFreq(site), GetNstate());
            // state =
            // Random::DrawFromDiscreteDistribution(GetRootFreq(site),GetNstate());
        } else {
            RootPosteriorDraw(site);
        }
    }
    for (const Link *link = from->Next(); link != from; link = link->Next()) {
        GetState(link->Out()->GetNode(), site) =
            GetSubMatrix(link->GetBranch()->GetIndex(), site)
                .DrawFiniteTime(state,
                                GetBranchLength(link->GetBranch()->GetIndex()) * GetSiteRate(site));
        PriorSample(link->Out(), site, rootprior);
    }
}

void PhyloProcess::ResampleState() {
    for (int i = 0; i < GetNsite(); i++) {
        ResampleState(i);
    }
}

void PhyloProcess::ResampleState(int site) {
    Pruning(GetRoot(), site);
    PruningAncestral(GetRoot(), site);
    // give information about fixed states at the tips to polyprocess
}

double PhyloProcess::Move(double fraction) {
    DrawSites(fraction);
    ResampleSub();
    // restoring full resampling mode
    DrawSites(1.0);
    return 1.0;
}

void PhyloProcess::DrawSites(double fraction) {
    for (int i = 0; i < GetNsite(); i++) {
        sitearray[i] = (Random::Uniform() < fraction);
    }
}

void PhyloProcess::ResampleSub() {
    pruningchrono.Start();
#if DEBUG > 1
    MeasureTime timer;
#endif

    for (int i = 0; i < GetNsite(); i++) {
        if (sitearray[i] != 0) {
            ResampleState(i);
        }
    }
#if DEBUG > 1
    timer.print<2>("ResampleSub - state. ");
#endif
    pruningchrono.Stop();

    resamplechrono.Start();
    for (int i = 0; i < GetNsite(); i++) {
        if (sitearray[i] != 0) {
            ResampleSub(GetRoot(), i);
        }
    }
    resamplechrono.Stop();
}

void PhyloProcess::ResampleSub(int site) {
    ResampleState(site);
    ResampleSub(GetRoot(), site);
}

void PhyloProcess::ResampleSub(const Link *from, int site) {
    delete pathmap[from->GetNode()][site];

    if (from->isRoot()) {
        pathmap[from->GetNode()][site] = SampleRootPath(GetState(from->GetNode(), site));
    } else {
        pathmap[from->GetNode()][site] =
            SamplePath(GetState(from->Out()->GetNode(), site), GetState(from->GetNode(), site),
                       GetBranchLength(from->GetBranch()->GetIndex()), GetSiteRate(site),
                       GetSubMatrix(from->GetBranch()->GetIndex(), site));
    }

    for (const Link *link = from->Next(); link != from; link = link->Next()) {
        ResampleSub(link->Out(), site);
    }
}

void PhyloProcess::PostPredSample(string name, bool rootprior) {
    for (int i = 0; i < GetNsite(); i++) {
        PostPredSample(i, rootprior);
    }
    SequenceAlignment tmpdata(*GetData());
    GetLeafData(&tmpdata);
    ofstream os(name.c_str());
    tmpdata.ToStream(os);
    os.close();
}

SequenceAlignment* PhyloProcess::PostPredSample(bool rootprior) {
    for (int i = 0; i < GetNsite(); i++) {
        PostPredSample(i, rootprior);
    }
    SequenceAlignment* tmpdata = new SequenceAlignment(*GetData());
    GetLeafData(tmpdata);
    return tmpdata;
}


void PhyloProcess::PostPredSample(int site, bool rootprior) {
    if (!rootprior) {
        Pruning(GetRoot(), site);
    }
    PriorSample(GetRoot(), site, rootprior);
}

void PhyloProcess::GetLeafData(SequenceAlignment *data) { RecursiveGetLeafData(GetRoot(), data); }

void PhyloProcess::RecursiveGetLeafData(const Link *from, SequenceAlignment *data) {
    if (from->isLeaf()) {
        for (int site = 0; site < GetNsite(); site++) {
            int state = GetState(from->GetNode(), site);
            int obsstate = GetData(from->GetNode()->GetIndex(), site);
            if (obsstate != unknown) {
                data->SetState(from->GetNode()->GetIndex(), site, state);
            } else {
                data->SetState(from->GetNode()->GetIndex(), site, unknown);
            }
        }
    }
    for (const Link *link = from->Next(); link != from; link = link->Next()) {
        RecursiveGetLeafData(link->Out(), data);
    }
}

BranchSitePath *PhyloProcess::SampleRootPath(int rootstate) {
    BranchSitePath *path = new BranchSitePath(rootstate);
    return path;
}

BranchSitePath *PhyloProcess::SamplePath(int stateup, int statedown, double time, double rate,
                                         const SubMatrix &matrix) {
    BranchSitePath *path = ResampleAcceptReject(1000, stateup, statedown, rate, time, matrix);
    if (!path) {
        path = ResampleUniformized(stateup, statedown, rate, time, matrix);
    }
    return path;
}

BranchSitePath *PhyloProcess::ResampleAcceptReject(int maxtrial, int stateup, int statedown,
                                                   double rate, double totaltime,
                                                   const SubMatrix &matrix) {
    int ntrial = 0;
    BranchSitePath *path = 0;

    if (rate * totaltime < 1e-10) {
        // if (rate * totaltime == 0)	{
        if (stateup != statedown) {
            cerr << "error in MatrixSubstitutionProcess::ResampleAcceptReject: "
                    "stateup != statedown, efflength == 0\n";
            exit(1);
        }
        delete path;
        path = new BranchSitePath();
        ntrial++;
        path->Reset(stateup);
    } else {
        do {
            delete path;
            path = new BranchSitePath();
            ntrial++;
            path->Reset(stateup);
            double t = 0;
            int state = stateup;

            if (state != statedown) {
                // draw waiting time conditional on at least one substitution
                double q = -rate * matrix(state, state);
                double u = -log(1 - Random::Uniform() * (1 - exp(-q * totaltime))) / q;

                t += u;
                int newstate = matrix.DrawOneStep(state);
                path->Append(newstate, u / totaltime);
                state = newstate;
            }
            while (t < totaltime) {
                // draw waiting time
                double q = -rate * matrix(state, state);
                double u = -log(1 - Random::Uniform()) / q;
                if (std::isnan(u)) {
                    cerr << "in MatrixSubstitutionProcess:: drawing exponential number: "
                            "nan\n";
                    cerr << rate << '\t' << q << '\n';
                    exit(1);
                }

                if (std::isinf(u)) {
                    cerr << "in MatrixSubstitutionProcess:: drawing exponential number: "
                            "inf\n";
                    cerr << rate << '\t' << q << '\n';
                    cerr << totaltime << '\n';
                    cerr << state << '\t' << stateup << '\n';
                    exit(1);
                }

                t += u;
                if (t < totaltime) {
                    int newstate = matrix.DrawOneStep(state);
                    path->Append(newstate, u / totaltime);
                    state = newstate;
                } else {
                    t -= u;
                    u = totaltime - t;
                    path->Last()->SetRelativeTime(u / totaltime);
                    t = totaltime;
                }
            }
        } while ((ntrial < maxtrial) && (path->Last()->GetState() != statedown));
    }

    // if endstate does not match state at the corresponding end of the branch
    // just force it to match
    // however, this is really dirty !
    // normally, in that case, one should give up with accept-reject
    // and use a uniformized method instead (but not yet adapted to the present
    // code, see below)
    if (path->Last()->GetState() != statedown) {
        // fossil
        // path->Last()->SetState(statedown);
        delete path;
        path = 0;
    }

    return path;
}

BranchSitePath *PhyloProcess::ResampleUniformized(int stateup, int statedown, double rate,
                                                  double totaltime, const SubMatrix &matrix) {
    double length = rate * totaltime;
    int m = matrix.DrawUniformizedSubstitutionNumber(stateup, statedown, length);

    vector<double> y(m + 1);
    for (int r = 0; r < m; r++) {
        y[r] = Random::Uniform();
    }
    y[m] = 1;
    sort(y.begin(), y.end());

    int state = stateup;

    BranchSitePath *path = new BranchSitePath();
    path->Reset(stateup);

    double t = y[0];
    for (int r = 0; r < m; r++) {
        int k = (r == m - 1) ? statedown
                             : matrix.DrawUniformizedTransition(state, statedown, m - r - 1);
        if (k != state) {
            path->Append(k, t);
            t = 0;
        }
        state = k;
        t += y[r + 1] - y[r];
    }
    path->Last()->SetRelativeTime(t);
    return path;
}

void PhyloProcess::AddPathSuffStat(PathSuffStat &suffstat) const {
    RecursiveAddPathSuffStat(GetRoot(), suffstat);
}

void PhyloProcess::RecursiveAddPathSuffStat(const Link *from, PathSuffStat &suffstat) const {
    LocalAddPathSuffStat(from, suffstat);
    for (const Link *link = from->Next(); link != from; link = link->Next()) {
        RecursiveAddPathSuffStat(link->Out(), suffstat);
    }
}

void PhyloProcess::LocalAddPathSuffStat(const Link *from, PathSuffStat &suffstat) const {
    int nodeindex = from->GetNode()->GetIndex();
    for (int i = 0; i < GetNsite(); i++) {
        if (missingmap[nodeindex][i] == 2) {
            suffstat.IncrementRootCount(GetState(from->GetNode(), i));
        } else if (missingmap[nodeindex][i] == 1) {
            if (from->isRoot()) {
                cerr << "error in missing map\n";
                exit(1);
            }
            pathmap[from->GetNode()][i]->AddPathSuffStat(
                suffstat, GetBranchLength(from->GetBranch()->GetIndex()) * GetSiteRate(i));
        }
    }
}

void PhyloProcess::AddPathSuffStat(BidimArray<PathSuffStat> &suffstatarray,
                                   const BranchAllocationSystem &branchalloc) const {
    RecursiveAddPathSuffStat(GetRoot(), suffstatarray, branchalloc);
}

void PhyloProcess::RecursiveAddPathSuffStat(const Link *from,
                                            BidimArray<PathSuffStat> &suffstatarray,
                                            const BranchAllocationSystem &branchalloc) const {
    if (from->isRoot()) {
        LocalAddPathSuffStat(from, suffstatarray, 0);
    } else {
        LocalAddPathSuffStat(from, suffstatarray,
                             branchalloc.GetBranchAlloc(from->GetBranch()->GetIndex()));
    }
    for (const Link *link = from->Next(); link != from; link = link->Next()) {
        RecursiveAddPathSuffStat(link->Out(), suffstatarray, branchalloc);
    }
}

void PhyloProcess::LocalAddPathSuffStat(const Link *from, BidimArray<PathSuffStat> &suffstatarray,
                                        int cond) const {
    int nodeindex = from->GetNode()->GetIndex();
    for (int i = 0; i < GetNsite(); i++) {
        if (missingmap[nodeindex][i] == 2) {
            suffstatarray(cond, i).IncrementRootCount(GetState(from->GetNode(), i));
        } else if (missingmap[nodeindex][i] == 1) {
            if (from->isRoot()) {
                cerr << "error in missing map\n";
                exit(1);
            }
            pathmap[from->GetNode()][i]->AddPathSuffStat(
                suffstatarray(cond, i),
                GetBranchLength(from->GetBranch()->GetIndex()) * GetSiteRate(i));
        }
    }
}

void PhyloProcess::AddPathSuffStat(Array<PathSuffStat> &suffstatarray) const {
    RecursiveAddPathSuffStat(GetRoot(), suffstatarray);
}

void PhyloProcess::RecursiveAddPathSuffStat(const Link *from,
                                            Array<PathSuffStat> &suffstatarray) const {
    LocalAddPathSuffStat(from, suffstatarray);
    for (const Link *link = from->Next(); link != from; link = link->Next()) {
        RecursiveAddPathSuffStat(link->Out(), suffstatarray);
    }
}

void PhyloProcess::LocalAddPathSuffStat(const Link *from,
                                        Array<PathSuffStat> &suffstatarray) const {
    int nodeindex = from->GetNode()->GetIndex();
    for (int i = 0; i < GetNsite(); i++) {
        if (missingmap[nodeindex][i] == 2) {
            suffstatarray[i].IncrementRootCount(GetState(from->GetNode(), i));
        } else if (missingmap[nodeindex][i] == 1) {
            if (from->isRoot()) {
                cerr << "error in missing map\n";
                exit(1);
            }
            pathmap[from->GetNode()][i]->AddPathSuffStat(
                suffstatarray[i], GetBranchLength(from->GetBranch()->GetIndex()) * GetSiteRate(i));
        }
    }
}

void PhyloProcess::AddPathSuffStat(NodeArray<PathSuffStat> &suffstatarray) const {
    RecursiveAddPathSuffStat(GetRoot(), suffstatarray);
}

void PhyloProcess::RecursiveAddPathSuffStat(const Link *from,
                                            NodeArray<PathSuffStat> &suffstatarray) const {
    LocalAddPathSuffStat(from, suffstatarray);
    for (const Link *link = from->Next(); link != from; link = link->Next()) {
        RecursiveAddPathSuffStat(link->Out(), suffstatarray);
    }
}

void PhyloProcess::LocalAddPathSuffStat(const Link *from,
                                        NodeArray<PathSuffStat> &suffstatarray) const {
    int nodeindex = from->GetNode()->GetIndex();
    for (int i = 0; i < GetNsite(); i++) {
        if (missingmap[nodeindex][i] == 2) {
            suffstatarray[nodeindex].IncrementRootCount(GetState(from->GetNode(), i));
        } else if (missingmap[nodeindex][i] == 1) {
            if (from->isRoot()) {
                cerr << "error in missing map\n";
                exit(1);
            }
            pathmap[from->GetNode()][i]->AddPathSuffStat(
                suffstatarray[nodeindex],
                GetBranchLength(from->GetBranch()->GetIndex()) * GetSiteRate(i));
        }
    }
}

void PhyloProcess::AddLengthSuffStat(
    BranchArray<PoissonSuffStat> &branchlengthpathsuffstatarray) const {
    RecursiveAddLengthSuffStat(GetRoot(), branchlengthpathsuffstatarray);
}

void PhyloProcess::RecursiveAddLengthSuffStat(
    const Link *from, BranchArray<PoissonSuffStat> &branchlengthpathsuffstatarray) const {
    if (!from->isRoot()) {
        LocalAddLengthSuffStat(from, branchlengthpathsuffstatarray[from->GetBranch()->GetIndex()]);
    }
    for (const Link *link = from->Next(); link != from; link = link->Next()) {
        RecursiveAddLengthSuffStat(link->Out(), branchlengthpathsuffstatarray);
    }
}

void PhyloProcess::LocalAddLengthSuffStat(const Link *link, PoissonSuffStat &suffstat) const {
    int nodeindex = link->GetNode()->GetIndex();
    for (int i = 0; i < GetNsite(); i++) {
        if (missingmap[nodeindex][i] == 1) {
            pathmap[link->GetNode()][i]->AddLengthSuffStat(
                suffstat, GetSiteRate(i), GetSubMatrix(link->GetBranch()->GetIndex(), i));
        }
    }
}

void PhyloProcess::AddRateSuffStat(Array<PoissonSuffStat> &siteratepathsuffstatarray) const {
    RecursiveAddRateSuffStat(GetRoot(), siteratepathsuffstatarray);
}

void PhyloProcess::RecursiveAddRateSuffStat(
    const Link *from, Array<PoissonSuffStat> &siteratepathsuffstatarray) const {
    if (!from->isRoot()) {
        LocalAddRateSuffStat(from, siteratepathsuffstatarray);
    }
    for (const Link *link = from->Next(); link != from; link = link->Next()) {
        RecursiveAddRateSuffStat(link->Out(), siteratepathsuffstatarray);
    }
}

void PhyloProcess::LocalAddRateSuffStat(const Link *link,
                                        Array<PoissonSuffStat> &siteratepathsuffstatarray) const {
    int nodeindex = link->GetNode()->GetIndex();
    double length = GetBranchLength(link->GetBranch()->GetIndex());
    for (int i = 0; i < GetNsite(); i++) {
        if (missingmap[nodeindex][i] == 1) {
            pathmap[link->GetNode()][i]->AddLengthSuffStat(
                siteratepathsuffstatarray[i], length,
                GetSubMatrix(link->GetBranch()->GetIndex(), i));
        }
    }
}
