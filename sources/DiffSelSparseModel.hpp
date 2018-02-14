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
#include "MultiGammaSuffStat.hpp"
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

    int codonmodel;

    // 0: free wo shrinkage
    // 1: free with shrinkage
    // 2: shared across genes
    // 3: fixed

    int blmode;
    int nucmode;
    int fitnesshypermode;

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

    // shiftprob (across conditions):
    // either Beta(shiftprobhypermean,shiftprobhyperinvconc), estimated across genes
    // or mixture 1-pi * 0 + pi * Beta: and this, for each condition separately
    vector<double> pi;
    vector<double> shiftprobhypermean;
    vector<double> shiftprobhyperinvconc;
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
    MultiGammaSuffStat hyperfitnesssuffstat;

  public:

    DiffSelSparseModel(const std::string& datafile, const std::string& treefile, int inNcond, int inNlevel, int incodonmodel) : hyperfitnesssuffstat(Naa) {

        codonmodel = incodonmodel;

        blmode = 0;
        nucmode = 0;
        fitnesshypermode = 3;

        Ncond = inNcond;
        Nlevel = inNlevel;

        ReadFiles(datafile, treefile);

        // specifies which condition for which branch
        branchalloc = new BranchAllocationSystem(*tree,Ncond);
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

        fitnessshape = 4.0;
        fitnesscenter.assign(Naa,1.0/Naa);
        fitness = new BidimIIDMultiGamma(Ncond,Nsite,Naa,fitnessshape,fitnesscenter);

        pi.assign(Ncond-1,0.1);
        shiftprobhypermean.assign(Ncond-1,0.5);
        shiftprobhyperinvconc.assign(Ncond-1,0.5);
        shiftprob.assign(Ncond-1,0.1);

        toggle = new BidimIIDMultiBernoulli(Ncond-1,Nsite,Naa,shiftprob);

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
        phyloprocess->Unfold();

        // create suffstat arrays
        suffstatarray = new PathSuffStatBidimArray(Ncond,Nsite);
    }

    void SetBLMode(int in)   {
        blmode = in;
    }

    void SetNucMode(int in) {
        nucmode = in;
    }

    void SetFitnessHyperMode(int in)    {
        fitnesshypermode = in;
    }

    // ------------------
    // Update system
    // ------------------

    //! \brief set branch lengths to a new value
    //! 
    //! Used in a multigene context.
    void SetBranchLengths(const BranchSelector<double>& inbranchlength)    {
        branchlength->Copy(inbranchlength);
    }

    void SetShiftProbHyperParameters(const vector<double>& inpi, const vector<double>& inshiftprobhypermean, const vector<double>& inshiftprobhyperinvconc)    {
        pi = inpi;
        shiftprobhypermean = inshiftprobhypermean;
        shiftprobhyperinvconc = inshiftprobhyperinvconc;
    }

    const vector<double>& GetShiftProbVector() const {
        return shiftprob;
    }


    void GetFitnessArray(int k, double* array) const  {
        int j = 0;
        for (int i=0; i<GetNsite(); i++)    {
            for (int a=0; a<Naa; a++)   {
                array[j++] = fitness->GetVal(k,i)[a];
            }
        }
    }

    void GetShiftToggleArray(int k, int* array) const  {
        int j = 0;
        for (int i=0; i<GetNsite(); i++)    {
            for (int a=0; a<Naa; a++)   {
                array[j++] = toggle->GetVal(k-1,i)[a];
            }
        }
    }

    void Update() override {
        branchlength->SetScale(lambda);
		fitness->SetShape(fitnessshape);
        UpdateAll();
        ResampleSub(1.0);
    }

    void NoUpdate() {}

    void CorruptMatrices()  {
        CorruptNucMatrix();
        condsubmatrixarray->Corrupt();
    }

    void CorruptNucMatrix() {
        nucmatrix->CopyStationary(nucstat);
        nucmatrix->CorruptMatrix();
    }

    void UpdateAll() {
        fitnessprofile->Update();
        CorruptMatrices();
    }

    void UpdateSite(int i) {
        fitnessprofile->UpdateColumn(i);
        condsubmatrixarray->CorruptColumn(i);
    }

    // ---------------
    // log priors
    // ---------------

    double GetLogPrior() const {
        double total = 0;
        if (blmode < 2) {
            total += BranchLengthsHyperLogPrior();
            total += BranchLengthsLogPrior();
        }
        if (nucmode < 2)    {
            total += NucRatesLogPrior();
        }
        if (fitnesshypermode < 2)   {
            total += FitnessHyperLogPrior();
        }
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
        double total = 0;
        for (int k=1; k<Ncond; k++) {
            double alpha = shiftprobhypermean[k-1] / shiftprobhyperinvconc[k-1];
            double beta = (1-shiftprobhypermean[k-1]) / shiftprobhyperinvconc[k-1];
            if (shiftprob[k-1] != 0)    {
                total += log(pi[k-1]) + Random::logBetaDensity(shiftprob[k-1],alpha,beta);
            }
            else    {
                total += log(1 - pi[k-1]);
            }
        }
        return total;
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

    //! \brief const access to array of length-pathsuffstats across branches
    //!
    //! Useful for resampling branch lengths conditional on the current substitution mapping
    const PoissonSuffStatBranchArray* GetLengthPathSuffStatArray() const {
        return lengthpathsuffstatarray;
    }

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

	double BranchLengthsHyperSuffStatLogProb()	const {
		return hyperlengthsuffstat.GetLogProb(1.0,lambda);
	}

	double FitnessHyperSuffStatLogProb()	const {
		return hyperfitnesssuffstat.GetLogProb(fitnessshape,fitnesscenter);
	}

    int GetNshift(int cond) const {
        if (! cond) {
            cerr << "error: GetNshift called on baseline\n";
            exit(1);
        }
        return toggle->GetRowEventNumber(cond-1);
    }

    int GetNshift(int cond, int site) const {
        if (! cond) {
            cerr << "error: GetNshift called on baseline\n";
            exit(1);
        }
        return toggle->GetEventNumber(cond-1,site);
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

    double FitnessHyperLogProb() const  {
        return FitnessHyperLogPrior() + FitnessHyperSuffStatLogProb();
    }

    // ---------------
    // Moves
    // ---------------

    //! \brief complete MCMC move schedule
	double Move() override {
        ResampleSub(1.0);
        MoveParameters(3,20);
        return 1.0;
	}

    void MoveParameters(int nrep0, int nrep) {

        for (int rep0 = 0; rep0 < nrep0; rep0++) {

            if (blmode < 2)    {
                CollectLengthSuffStat();
                ResampleBranchLengths();
                MoveBranchLengthsHyperParameter();
            }

            CollectPathSuffStat();
            UpdateAll();

            for (int rep = 0; rep < nrep; rep++) {
                MoveBaselineFitness();
                CompMoveFitness();
                MoveFitnessShifts();
                MoveShiftToggles();
                if (fitnesshypermode < 2)   {
                    MoveFitnessHyperParameters();
                }
            }

            if (nucmode < 2)    {
                MoveNucRates();
            }
        }

        UpdateAll();
    }

    void ResampleSub(double frac)   {
        CorruptMatrices();
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
        MoveBaselineFitness(1.0, 3, 10);
        MoveBaselineFitness(1.0, 10, 10);
        MoveBaselineFitness(1.0, 20, 10);
        MoveBaselineFitness(0.3, 20, 10);
    }

    void CompMoveFitness()  {
        CompMoveFitness(1.0,10);
    }

    double CompMoveFitness(double tuning, int nrep) {

        double nacc = 0;
        double ntot = 0;

        for (int rep = 0; rep < nrep; rep++) {
            for (int i = 0; i < Nsite; i++) {

                double deltalogprob = 0;

                for (int k=0; k<Ncond; k++) {
                    for (int a=0; a<Naa; a++)   {
                        if ((!k) || ((*toggle)(k-1,i)[a]))   {
                            double alpha = fitnessshape*fitnesscenter[a];
                            deltalogprob -= - Random::logGamma(alpha) + (alpha-1)*log((*fitness)(k,i)[a]) - (*fitness)(k,i)[a];
                        }
                    }
                }

                double m = tuning*(Random::Uniform() - 0.5);
                double e = exp(m);

                int n = 0;
                for (int k=0; k<Ncond; k++) {
                    for (int a=0; a<Naa; a++)   {
                        if ((!k) || ((*toggle)(k-1,i)[a]))   {
                            (*fitness)(k,i)[a] *= e;
                            n++;
                        }
                    }
                }

                double loghastings = n * m;

                for (int k=0; k<Ncond; k++) {
                    for (int a=0; a<Naa; a++)   {
                        if ((!k) || ((*toggle)(k-1,i)[a]))   {
                            double alpha = fitnessshape*fitnesscenter[a];
                            deltalogprob += - Random::logGamma(alpha) + (alpha-1)*log((*fitness)(k,i)[a]) - (*fitness)(k,i)[a];
                        }
                    }
                }

                deltalogprob += loghastings;

                int accepted = (log(Random::Uniform()) < deltalogprob);
                if (accepted) {
                    nacc++;
                } else {
                    for (int k=0; k<Ncond; k++) {
                        for (int a=0; a<Naa; a++)   {
                            if ((!k) || ((*toggle)(k-1,i)[a]))   {
                                (*fitness)(k,i)[a] /= e;
                            }
                        }
                    }
                }
                ntot++;
            }
        }
        return nacc / ntot;
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
            MoveFitnessShifts(k,1,10);
            MoveFitnessShifts(k,0.3,10);
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

    void MoveFitnessHyperParameters() {
        // collect suff stats across all active fitness parameters
        hyperfitnesssuffstat.Clear();
        hyperfitnesssuffstat.AddSuffStat(*fitness,*toggle);

        ScalingMove(fitnessshape,1.0,100,&DiffSelSparseModel::FitnessHyperLogProb,&DiffSelSparseModel::NoUpdate,this);
        ScalingMove(fitnessshape,0.3,100,&DiffSelSparseModel::FitnessHyperLogProb,&DiffSelSparseModel::NoUpdate,this);
        ScalingMove(fitnessshape,0.1,100,&DiffSelSparseModel::FitnessHyperLogProb,&DiffSelSparseModel::NoUpdate,this);

        ProfileMove(fitnesscenter,0.3,1,100,&DiffSelSparseModel::FitnessHyperLogProb,&DiffSelSparseModel::NoUpdate,this);
        ProfileMove(fitnesscenter,0.1,1,100,&DiffSelSparseModel::FitnessHyperLogProb,&DiffSelSparseModel::NoUpdate,this);
        ProfileMove(fitnesscenter,0.1,3,100,&DiffSelSparseModel::FitnessHyperLogProb,&DiffSelSparseModel::NoUpdate,this);

		fitness->SetShape(fitnessshape);

    }

    void ResampleShiftProb()    {

        for (int k=1; k<Ncond; k++) {

            double alpha = shiftprobhypermean[k-1] / shiftprobhyperinvconc[k-1];
            double beta = (1-shiftprobhypermean[k-1]) / shiftprobhyperinvconc[k-1];

            int nshift = 0;
            for (int i=0; i<Nsite; i++) {
                for (int a=0; a<Naa; a++)   {
                    nshift += (*toggle)(k-1,i)[a];
                }
            }
            int nn = Nsite*Naa;

            if (nshift) {
                shiftprob[k-1] = Random::BetaSample(alpha + nshift, beta + nn - nshift);
            }
            else    {
                double logp0 = log(1-pi[k-1]);
                double logp1 = log(pi[k-1]) + Random::logGamma(alpha+beta) + Random::logGamma(beta + nn) - Random::logGamma(beta) - Random::logGamma(alpha+beta+nn);
                double max = (logp0 > logp1) ? logp0 : logp1;
                double p0 = exp(logp0-max);
                double p1 = exp(logp1-max);
                double tot = p0+p1;
                p0/=tot;
                p1/=tot;
                if (Random::Uniform() < p0) {
                    shiftprob[k-1] = 0;
                }
                else    {
                    shiftprob[k-1] = Random::BetaSample(alpha + nshift, beta + nn - nshift);
                }
            }
        }
    }

    void MoveShiftToggles() {
        for (int k=1; k<Ncond; k++) {
            MoveShiftToggles(k,10);
        }
    }

    double ToggleMarginalLogPrior(int nn, int nshift, double pi, double alpha, double beta) const {
        double logp1 = log(pi) + Random::logGamma(alpha+beta) - Random::logGamma(alpha) - Random::logGamma(beta) + Random::logGamma(alpha+nshift) + Random::logGamma(beta + nn - nshift) - Random::logGamma(alpha + beta + nn);
        if (nshift) {
            return logp1;
        }
        double logp0 = log(1-pi);
        double max = (logp0 > logp1) ? logp0 : logp1;
        double tot = exp(logp0-max) + exp(logp1-max);
        double ret = log(tot) + max;
        return ret;
    }

    double MoveShiftToggles(int k, int nrep)  {

        int nshift = 0;
        for (int i=0; i<Nsite; i++) {
            for (int a=0; a<Naa; a++)   {
                nshift += (*toggle)(k-1,i)[a];
            }
        }
        int nn = Nsite*Naa;

        double alpha = shiftprobhypermean[k-1] / shiftprobhyperinvconc[k-1];
        double beta = (1-shiftprobhypermean[k-1]) / shiftprobhyperinvconc[k-1];
        double pp = pi[k-1];

        double ntot = 0;
        double nacc = 0;
        for (int rep = 0; rep < nrep; rep++) {
            for (int i = 0; i < Nsite; i++) {
                int a = (int) (Naa * Random::Uniform());

                if (!(*toggle)(k-1,i)[a])    {

                    double deltalogprob = -ToggleMarginalLogPrior(nn,nshift,pp,alpha,beta) - SiteSuffStatLogProb(i);
                    (*toggle)(k-1,i)[a] = 1;
                    (*fitness)(k,i)[a] = Random::sGamma(fitnessshape * fitnesscenter[a]);
                    // (*fitness)(k,i)[a] = Random::Gamma(fitnessshape, fitnessshape / fitnesscenter[a]);
                    UpdateSite(i);
                    deltalogprob += ToggleMarginalLogPrior(nn,nshift+1,pp,alpha,beta) + SiteSuffStatLogProb(i);
                    /*
                    double c1 = ToggleMarginalLogPrior(nn,nshift+1,pp,alpha,beta) - ToggleMarginalLogPrior(nn,nshift,pp,alpha,beta);
                    double c2 = log(alpha + nshift) - log(beta + nn - nshift - 1);
                    cerr << nshift << c2 - c1 << '\n';
                    */
                    // deltalogprob += log(alpha + nshift) - log(beta + nn - nshift - 1);

                    int accepted = (log(Random::Uniform()) < deltalogprob);
                    if (accepted) {
                        nacc++;
                        nshift++;
                    } else {
                        (*toggle)(k-1,i)[a] = 0;
                        UpdateSite(i);
                    }
                    ntot++;
                }
                else    {
                    double deltalogprob = -ToggleMarginalLogPrior(nn,nshift,pp,alpha,beta) - SiteSuffStatLogProb(i);
                    (*toggle)(k-1,i)[a] = 0;
                    UpdateSite(i);
                    deltalogprob += ToggleMarginalLogPrior(nn,nshift-1,pp,alpha,beta) + SiteSuffStatLogProb(i);
                    // deltalogprob += log(beta + nn - nshift) + log(alpha + nshift - 1);

                    int accepted = (log(Random::Uniform()) < deltalogprob);
                    if (accepted) {
                        nacc++;
                        nshift--;
                    } else {
                        (*toggle)(k-1,i)[a] = 1;
                        UpdateSite(i);
                    }
                    ntot++;
                }
            }
        }
        shiftprob[k-1] = Random::BetaSample(alpha + nshift, beta + nn - nshift);
        return nacc / ntot;
    }

    //-------------------
    // Accessors
    // ------------------

	CodonStateSpace* GetCodonStateSpace() const {
		return (CodonStateSpace*) codondata->GetStateSpace();
	}

    int GetNsite() const { return Nsite; }
    int GetNcond() const { return Ncond; }

    //-------------------
    // Traces and monitors
    // ------------------

    void TraceHeader(std::ostream& os) const override {
        os << "#logprior\tlnL\tlength\t";
        os << "meanvar0\t";
        os << "shape\t";
        os << "center\t";
        for (int k = 1; k < Ncond; k++) {
            os << "prob" << k << '\t';
        }
        os << "statent\t";
        os << "rrent\n";
    }

    void Trace(ostream& os) const override {
        os << GetLogPrior() << '\t';
        os << GetLogLikelihood() << '\t';
        os << branchlength->GetTotalLength() << '\t';
        os << fitness->GetMeanRelVar(0) << '\t';
        os << fitnessshape << '\t';
        os << Random::GetEntropy(fitnesscenter) << '\t';
        for (int k=1; k<Ncond; k++) {
            os << shiftprob[k-1] << '\t';
        }
        os << Random::GetEntropy(nucstat) << '\t';
        os << Random::GetEntropy(nucrelrate) << '\n';
    }

    void TraceToggle(int k, ostream& os) const {
        for (int i=0; i<GetNsite(); i++)    {
            for (int a=0; a<Naa; a++)   {
                os << toggle->GetVal(k-1,i)[a] << '\t';
            }
        }
        os << '\n';
    }

    void TraceFitness(int k, ostream& os) const {
        for (int i=0; i<GetNsite(); i++)    {
            for (int a=0; a<Naa; a++)   {
                os << fitness->GetVal(k,i)[a] << '\t';
            }
        }
        os << '\n';
    }

    void Monitor(ostream&) const override {}

    void FromStream(istream& is) override {
        if (blmode < 2) {
            is >> lambda;
            is >> *branchlength;
        }
        if (nucmode < 2)    {
            is >> nucrelrate;
            is >> nucstat;
        }
        if (fitnesshypermode < 2)   {
            is >> fitnessshape;
            is >> fitnesscenter;
        }
        is >> *fitness;
        is >> shiftprob;
        is >> *toggle;
    }

    void ToStream(ostream& os) const override {
        if (blmode < 2) {
            os << lambda << '\t';
            os << *branchlength << '\t';
        }
        if (nucmode < 2)    {
            os << nucrelrate << '\t';
            os << nucstat << '\t';
        }
        if (fitnesshypermode < 2)   {
            os << fitnessshape << '\t';
            os << fitnesscenter << '\t';
        }
        os << *fitness << '\t';
        os << shiftprob << '\t';
        os << *toggle << '\t';
    }

    //! return size of model, when put into an MPI buffer (in multigene context -- only omegatree)
    unsigned int GetMPISize() const {
        int size = 0;
        if (blmode < 2) {
            size++;
            size += branchlength->GetMPISize();
        }
        if (nucmode < 2)    {
            size += nucrelrate.size();
            size += nucstat.size();
        }
        if (fitnesshypermode < 2)   {
            size ++;
            size += fitnesscenter.size();
        }
        size += fitness->GetMPISize();
        size += shiftprob.size();
        size += toggle->GetMPISize();
        return size;
    }

    //! get array from MPI buffer
    void MPIGet(const MPIBuffer& is)    {
        if (blmode < 2) {
            is >> lambda;
            is >> *branchlength;
        }
        if (nucmode < 2)    {
            is >> nucrelrate;
            is >> nucstat;
        }
        if (fitnesshypermode < 2)   {
            is >> fitnessshape;
            is >> fitnesscenter;
        }
        is >> *fitness;
        is >> shiftprob;
        is >> *toggle;
    }

    //! write array into MPI buffer
    void MPIPut(MPIBuffer& os) const {
        if (blmode < 2) {
            os << lambda;
            os << *branchlength;
        }
        if (nucmode < 2)    {
            os << nucrelrate;
            os << nucstat;
        }
        if (fitnesshypermode < 2)   {
            os << fitnessshape;
            os << fitnesscenter;
        }
        os << *fitness;
        os << shiftprob;
        os << *toggle;
    }
};
