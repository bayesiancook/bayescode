/*Copyright or © or Copr. Centre National de la Recherche Scientifique (CNRS) (2017-06-14).
Contributors:
* Nicolas LARTILLOT - nicolas.lartillot@univ-lyon1.fr

This software is a computer program whose purpose is to detect convergent evolution using Bayesian
phylogenetic codon models.

This software is governed by the CeCILL-C license under French law and abiding by the rules of
distribution of free software. You can use, modify and/ or redistribute the software under the terms
of the CeCILL-C license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info".

As a counterpart to the access to the source code and rights to copy, modify and redistribute
granted by the license, users are provided only with a limited warranty and the software's author,
the holder of the economic rights, and the successive licensors have only limited liability.

In this respect, the user's attention is drawn to the risks associated with loading, using,
modifying and/or developing or reproducing the software by the user in light of its specific status
of free software, that may mean that it is complicated to manipulate, and that also therefore means
that it is reserved for developers and experienced professionals having in-depth computer knowledge.
Users are therefore encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or data to be ensured and,
more generally, to use and operate it in the same conditions as regards security.

The fact that you are presently reading this means that you have had knowledge of the CeCILL-C
license and that you accept its terms.*/


#include "CodonSequenceAlignment.hpp"
#include "GTRSubMatrix.hpp"
#include "PhyloProcess.hpp"
#include "ProbModel.hpp"
#include "Tree.hpp"
#include "IIDMultiGamma.hpp"
#include "IIDMultiBernoulli.hpp"
#include "DiffSelSparseFitnessArray.hpp"
#include "BranchAllocationSystem.hpp"
#include "AAMutSelCodonMatrixArray.hpp"
#include "SubMatrixSelector.hpp"
#include "IIDGamma.hpp"
#include "GammaSuffStat.hpp"
#include "IIDDirichlet.hpp"
#include "PathSuffStat.hpp"

class DiffSelSparseModel : public ProbModel {

    // -----
    // model selectors
    // -----

    int fixglob;
    int fixvar;
    int codonmodel;

    // -----
    // external parameters
    // -----

    Tree* tree;
    FileSequenceAlignment* data;
    const TaxonSet* taxonset;
    CodonSequenceAlignment* codondata;

    // number of sites
    int Nsite;
    int Ntaxa;
    int Nbranch;

    // number of diff sel categories
    int Ncond;

    // number of levels of the model
    // with 2 levels, structure of the model is as follows:
    // baseline (condition 0)
    // baseline  || fitness1 (for condition 1)
    // baseline || fitness1  || fitnessk  (for condition k=2..Ncond)
    int Nlevel;

    // which branch is under which condition
    BranchAllocationSystem* branchalloc;

    // auxiliary array
    // condalloc[k][l] is 1 iff condition l should recompute its fitness profiles whenever condition k has changed
    // has changed
    vector<vector<int> > condalloc;

    // -----
    //  model structure
    // -----

    // branch lengths iid expo (gamma of shape 1 and scale lambda)
    // where lambda is a hyperparameter
	double lambda;
	BranchIIDGamma* branchlength;

    // nucleotide exchange rates and equilibrium frequencies (stationary probabilities)
	std::vector<double> nucrelrate;
	std::vector<double> nucstat;
	GTRSubMatrix* nucmatrix;

    double fitnessshape;
    vector<double> fitnesscenter;
    BidimIIDMultiGamma* fitness;

    vector<double> shiftprob;
    BidimIIDMultiBernoulli* toggle;

    // fitness profiles (combinations of baseline and delta)
    // across conditions and across sites
    DiffSelSparseFitnessArray* fitnessprofile;

    // codon substitution matrices
    // across conditions and sites
    AAMutSelCodonMatrixArray* condsubmatrixarray;

    // branch- and site-substitution matrices (for phyloprocess)
    SubMatrixSelector* submatrixarray;
    // and for root (condition 0)
    RootSubMatrixSelector* rootsubmatrixarray;

    // phyloprocess
    PhyloProcess* phyloprocess;

    // suff stats

    // path suff stats across conditions and sites
    PathSuffStatBidimArray* suffstatarray;

    // Poisson suffstats for substitution histories, as a function of branch lengths
	PoissonSuffStatBranchArray* lengthpathsuffstatarray;

    // suff stats branch lengths, as a function of their hyper parameter lambda
    // (bl are iid gamma, of scale parameter lambda)
	GammaSuffStat hyperlengthsuffstat;

  public:

    DiffSelSparseModel(const std::string& datafile, const std::string& treefile, int inNcond, int inNlevel, int incodonmodel) {

        codonmodel = incodonmodel;
        Ncond = inNcond;
        Nlevel = inNlevel;

        ReadFiles(datafile, treefile);
        MakePattern();

        // specifies which condition for which branch
        branchalloc = new BranchAllocationSystem(*tree,Ncond);
        std::cerr << "-- conditions over branches ok\n";
    }

    // DiffSelSparseModel(const DiffSelSparseModel&) = delete;

    ~DiffSelSparseModel() {}

    void ReadFiles(string datafile, string treefile) {
        // nucleotide sequence alignment
        data = new FileSequenceAlignment(datafile);

        // translated into codon sequence alignment
        codondata = new CodonSequenceAlignment(data, true);

        Nsite = codondata->GetNsite();  // # columns
        Ntaxa = codondata->GetNtaxa();

        std::cerr << "-- Number of sites: " << Nsite << std::endl;

        taxonset = codondata->GetTaxonSet();

        // get tree from file (newick format)
        tree = new Tree(treefile);

        // check whether tree and data fits together
        tree->RegisterWith(taxonset);

        // traversal of the tree, so as to number links, branches and nodes
        // convention is: branches start at 1 (branch number 0 is the null branch behind the root)
        // nodes start at 0 (for the root), and nodes 1..Ntaxa are tip nodes (corresponding to taxa
        // in sequence alignment)
        tree->SetIndices();
        Nbranch = tree->GetNbranch();

        std::cerr << "-- Number of taxa : " << Ntaxa << '\n';
        std::cerr << "-- Number of branches : " << Nbranch << '\n';

        std::cerr << "-- Tree and data fit together\n";
    }

    void MakePattern()  {
        condalloc.assign(Ncond,vector<int>(Ncond,0));
        for (int l = 0; l < Ncond; l++) {
            condalloc[0][l] = 1;
        }
        if (Nlevel == 2) {
            for (int l = 1; l < Ncond; l++) {
                condalloc[1][l] = 1;
            }
        }
        for (int l = Nlevel; l < Ncond; l++) {
            condalloc[l][l] = 1;
        }
    }

    void Allocate() {

        // ----------
        // construction of the model
        // ----------

        // allocating data structures and sampling initial configuration

        // branch lengths
		lambda = 10;
		branchlength = new BranchIIDGamma(*tree,1.0,lambda);
		lengthpathsuffstatarray = new PoissonSuffStatBranchArray(*tree);

        // nucleotide matrix
		nucrelrate.assign(Nrr,0);
        Random::DirichletSample(nucrelrate,vector<double>(Nrr,1.0/Nrr),((double) Nrr));
		nucstat.assign(Nnuc,0);
        Random::DirichletSample(nucstat,vector<double>(Nnuc,1.0/Nnuc),((double) Nnuc));
		nucmatrix = new GTRSubMatrix(Nnuc,nucrelrate,nucstat,true);

        fitnessshape = 1.0;
        fitnesscenter.assign(Naa,1.0/Naa);
        fitness = new BidimIIDMultiGamma(Nsite,Ncond,Naa,fitnessshape,fitnesscenter);

        shiftprob.assign(Ncond-1,0.1);
        toggle = new BidimIIDMultiBernoulli(Nsite,Ncond-1,Naa,shiftprob);

        fitnessprofile = new DiffSelSparseFitnessArray(*fitness,*toggle,Nlevel);
        
        // codon matrices
        // per condition and per site
        condsubmatrixarray = new AAMutSelCodonMatrixArray(*fitnessprofile,*GetCodonStateSpace(),*nucmatrix);

        // sub matrices per branch and per site
        submatrixarray = new SubMatrixSelector(*condsubmatrixarray,*branchalloc);
        // sub matrices for root, across sites
        rootsubmatrixarray = new RootSubMatrixSelector(*condsubmatrixarray);

        // create phyloprocess
        phyloprocess = new PhyloProcess(tree, codondata, branchlength, 0, submatrixarray, rootsubmatrixarray);

        // create suffstat arrays
        suffstatarray = new PathSuffStatBidimArray(Ncond,Nsite);
    }

    void Unfold(bool sample)   {

        // unfold phyloprocess (allocate conditional likelihood vectors, etc)
        std::cerr << "-- unfolding\n";
        phyloprocess->Unfold();

        if (sample) {
            // stochastic mapping of substitution histories
            std::cerr << "-- mapping substitutions\n";
            phyloprocess->ResampleSub();
        }
    }

    // ------------------
    // Update system
    // ------------------

    void NoUpdate() {}

    void CorruptMatrices()  {
        CorruptNucMatrix();
        condsubmatrixarray->Corrupt();
    }

    void CorruptNucMatrix() {
        nucmatrix->CopyStationary(nucstat);
        nucmatrix->CorruptMatrix();
    }

    void Update() override {
        fitnessprofile->Update();
        CorruptMatrices();
        phyloprocess->GetLogLikelihood();
    }

    void UpdateAll() {
        fitnessprofile->Update();
        CorruptMatrices();
    }

    void UpdateSite(int i) {
        fitnessprofile->UpdateColumn(i);
        condsubmatrixarray->CorruptColumn(i);
    }

    void UpdateSiteCond(int i, int k)   {
        fitnessprofile->UpdateColumn(i,condalloc[k+1]);
        condsubmatrixarray->CorruptColumn(i,condalloc[k+1]);
    }

    // ---------------
    // log priors
    // ---------------

    double GetLogPrior() const {
        double total = 0;

        // branchlengths
        total += BranchLengthsHyperLogPrior();
        total += BranchLengthsLogPrior();

        // nuc rates
        total += NucRatesLogPrior();

        total += FitnessHyperLogPrior();
        total += FitnessLogPrior();

        total += ToggleHyperLogPrior();
        total += ToggleLogPrior();

        return total;
    }

	double BranchLengthsHyperLogPrior()	const {
        // exponential of mean 10
		return -lambda / 10;
	}

	double BranchLengthsLogPrior()	const {
		return branchlength->GetLogProb();
	}

    double NucRatesLogPrior() const {
        // uniform on relrates and nucstat
        double total = 0;
        total += Random::logGamma((double)Nnuc);
        total += Random::logGamma((double)Nrr);
        return total;
    }

    double FitnessHyperLogPrior() const {
        // uniform on center
        // exponential on shape
        return -fitnessshape;
    }

    double FitnessLogPrior() const  {
        return fitness->GetLogProb();
    }

    // for the moment, uniform prior over shift probs
    double ToggleHyperLogPrior() const  {
        return 0;
    }

    double ToggleLogPrior() const   {
        return toggle->GetLogProb();
    }

    double GetLogLikelihood() const { 
        return phyloprocess->GetLogLikelihood();
    }

    double GetLogProb() const {
        return GetLogPrior() + GetLogLikelihood();
    }

    // ---------------
    // collecting suff stats
    // ---------------

    // suffstats, per condition and per site
    // see SuffStat.hpp
    void CollectPathSuffStat() {
        suffstatarray->Clear();
        suffstatarray->AddSuffStat(*phyloprocess,*branchalloc);
    }

    void CollectLengthSuffStat()    {
		lengthpathsuffstatarray->Clear();
        lengthpathsuffstatarray->AddLengthPathSuffStat(*phyloprocess);
    }

    double SuffStatLogProb() const   {
        return suffstatarray->GetLogProb(*condsubmatrixarray);
    }

    double SiteSuffStatLogProb(int site) const   {
        return suffstatarray->GetLogProb(site,*condsubmatrixarray);
    }

    double SiteCondSuffStatLogProb(int site,int k)  {
        return suffstatarray->GetLogProb(site,condalloc[k+1],*condsubmatrixarray);
    }

	double BranchLengthsHyperSuffStatLogProb()	const {
		return hyperlengthsuffstat.GetLogProb(1.0,lambda);
	}

    int GetNshift(int cond) const {
        if (! cond) {
            cerr << "error: GetNshift called on baseline\n";
            exit(1);
        }
        return toggle->GetRowEventNumber(cond);
    }

    int GetNshift(int cond, int site) const {
        if (! cond) {
            cerr << "error: GetNshift called on baseline\n";
            exit(1);
        }
        return toggle->GetEventNumber(cond,site);
    }

    // ---------------
    // log probs for MH moves
    // ---------------

    double BranchLengthsHyperLogProb() const {
        return BranchLengthsHyperLogPrior() + BranchLengthsHyperSuffStatLogProb();
    }

    double NucRatesLogProb() const {
        return NucRatesLogPrior() + SuffStatLogProb();
    }

    // ---------------
    // Moves
    // ---------------

    // move cycle schedule
    // does not yet implement any monitoring (success rates, time spent, etc)
    double Move() override {

        int nrep0 = 3;
        int nrep = 20;

        for (int rep0 = 0; rep0 < nrep0; rep0++) {

            CollectLengthSuffStat();

            ResampleBranchLengths();
            MoveBranchLengthsHyperParameter();

            CollectPathSuffStat();

            UpdateAll();

            for (int rep = 0; rep < nrep; rep++) {
                MoveBaselineFitness();
                MoveFitnessShifts();
                MoveShiftToggles();
                MoveFitnessHyperParameters();
                ResampleShiftProbs();
            }
            MoveNucRates();
        }

        UpdateAll();

        ResampleSub(1.0);
        return 1.0;
    }

    void MoveFitnessHyperParameters() {}
    void ResampleShiftProbs() {}
    void MoveShiftToggles() {}

    void ResampleSub(double frac)   {
		phyloprocess->Move(frac);
    }

	void ResampleBranchLengths()	{
        CollectLengthSuffStat();
		branchlength->GibbsResample(*lengthpathsuffstatarray);
	}

	void MoveBranchLengthsHyperParameter()	{

		hyperlengthsuffstat.Clear();
		hyperlengthsuffstat.AddSuffStat(*branchlength);
        ScalingMove(lambda,1.0,10,&DiffSelSparseModel::BranchLengthsHyperLogProb,&DiffSelSparseModel::NoUpdate,this);
        ScalingMove(lambda,0.3,10,&DiffSelSparseModel::BranchLengthsHyperLogProb,&DiffSelSparseModel::NoUpdate,this);
		branchlength->SetScale(lambda);
	}

	void MoveNucRates()	{

        CorruptMatrices();

        ProfileMove(nucrelrate,0.1,1,10,&DiffSelSparseModel::NucRatesLogProb,&DiffSelSparseModel::CorruptMatrices,this);
        ProfileMove(nucrelrate,0.03,3,10,&DiffSelSparseModel::NucRatesLogProb,&DiffSelSparseModel::CorruptMatrices,this);
        ProfileMove(nucrelrate,0.01,3,10,&DiffSelSparseModel::NucRatesLogProb,&DiffSelSparseModel::CorruptMatrices,this);

        ProfileMove(nucstat,0.1,1,10,&DiffSelSparseModel::NucRatesLogProb,&DiffSelSparseModel::CorruptMatrices,this);
        ProfileMove(nucstat,0.01,1,10,&DiffSelSparseModel::NucRatesLogProb,&DiffSelSparseModel::CorruptMatrices,this);

        CorruptMatrices();
	}

    void MoveBaselineFitness() {
        MoveBaselineFitness(0.15, 10, 1);
    }

    double MoveBaselineFitness(double tuning, int n, int nrep) {

        double nacc = 0;
        double ntot = 0;
        vector<double> bk(Naa,0);


        for (int rep = 0; rep < nrep; rep++) {
            for (int i = 0; i < Nsite; i++) {

                vector<double>& x = (*fitness)(0,i);

                bk = x;

                double deltalogprob = -fitness->GetLogProb(0,i) - SiteSuffStatLogProb(i);
                double loghastings = Random::PosRealVectorProposeMove(x, Naa, tuning, n);
                deltalogprob += loghastings;

                UpdateSite(i);

                deltalogprob += fitness->GetLogProb(0,i) + SiteSuffStatLogProb(i);

                int accepted = (log(Random::Uniform()) < deltalogprob);
                if (accepted) {
                    nacc++;
                } else {
                    x = bk;
                    UpdateSite(i);
                }
                ntot++;
            }
        }
        return nacc / ntot;
    }

    void MoveFitnessShifts()    {
        for (int k=1; k<Ncond; k++) {
            MoveFitnessShifts(k,1,3);
            MoveFitnessShifts(k,0.3,3);
        }
    }

    double MoveFitnessShifts(int k, double tuning, int nrep) {

        double nacc = 0;
        double ntot = 0;
        vector<double> bk(Naa,0);

        for (int rep = 0; rep < nrep; rep++) {
            for (int i = 0; i < Nsite; i++) {
                if (GetNshift(k,i))  {

                    vector<double>& x = (*fitness)(k,i);
                    const vector<int>& s = (*toggle)(k-1,i);

                    bk = x;

                    double deltalogprob = -fitness->GetLogProb(k,i,s) - SiteSuffStatLogProb(i);
                    double loghastings = Random::PosRealVectorProposeMove(x, Naa, tuning, s);
                    deltalogprob += loghastings;

                    UpdateSite(i);

                    deltalogprob += fitness->GetLogProb(k,i,s) + SiteSuffStatLogProb(i);

                    int accepted = (log(Random::Uniform()) < deltalogprob);
                    if (accepted) {
                        nacc++;
                    } else {
                        x = bk;
                        UpdateSite(i);
                    }
                    ntot++;
                }
            }
        }
        return nacc / ntot;
    }

    //-------------------
    // Accessors
    // ------------------

	CodonStateSpace* GetCodonStateSpace() const {
		return (CodonStateSpace*) codondata->GetStateSpace();
	}

    int GetNsite() { return Nsite; }
    int GetNcond() { return Ncond; }

    //-------------------
    // Traces and monitors
    // ------------------

    void TraceHeader(std::ostream& os) const override {
        os << "#logprior\tlnL\tlength\t";
        os << "globent\t";
        for (int k = 1; k < Ncond; k++) {
            os << "selvar" << k << '\t';
        }
        os << "statent\t";
        os << "rrent\n";
    }

    void Trace(ostream& os) const override {
        os << GetLogPrior() << '\t';
        os << GetLogLikelihood() << '\t';
        os << branchlength->GetTotalLength() << '\t';
        /*
        os << baseline->GetMeanEntropy() << '\t';
        for (int k = 0; k < Ncond-1; k++) {
            os << delta->GetMeanVar(k) << '\t';
        }
        */
        os << Random::GetEntropy(nucstat) << '\t';
        os << Random::GetEntropy(nucrelrate) << '\n';
    }

    void Monitor(ostream&) const override {}

    void FromStream(istream& is) override {
        is >> lambda;
        is >> *branchlength;
        is >> nucrelrate;
        is >> nucstat;
    }

    void ToStream(ostream& os) const override {
        os << lambda << '\n';
        os << *branchlength << '\n';
        os << nucrelrate << '\n';
        os << nucstat << '\n';
    }
};
