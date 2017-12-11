
#include "CodonSequenceAlignment.hpp"
#include "Tree.hpp"
#include "ProbModel.hpp"
#include "GTRSubMatrix.hpp"
#include "AAMutSelOmegaCodonSubMatrix.hpp"
#include "PhyloProcess.hpp"
#include "IIDGamma.hpp"
#include "IIDDirichlet.hpp"
#include "CodonSuffStat.hpp"
#include "ProbModel.hpp"
#include "Chrono.hpp"

class AAMutSelOmegaModel : public ProbModel {

	Tree* tree;
	FileSequenceAlignment* data;
	const TaxonSet* taxonset;
	CodonSequenceAlignment* codondata;

	int Nsite;
	int Ntaxa;
	int Nbranch;

	double lambda;
	BranchIIDGamma* branchlength;
	PoissonSuffStatBranchArray* lengthpathsuffstatarray;
	GammaSuffStat lambdasuffstat;

	std::vector<double> nucstat;
	std::vector<double> nucrelrate;
	GTRSubMatrix* nucmatrix;

    // of mean omegahypermean and inverse shape parameter omegahyperinvshape
    double omegahypermean;
    double omegahyperinvshape;
	double omega;
	OmegaSuffStat omegasuffstat;
	
    vector<double> aacenter;
    double aaconc;
    IIDDirichlet* aafitnessarray;
    IIDDirichlet* bkaafitnessarray;
	AAMutSelOmegaCodonSubMatrixArray* codonmatrixarray;
    DirichletSuffStat aahypersuffstat;

	PhyloProcess* phyloprocess;

	PathSuffStatArray* pathsuffstatarray;

    Chrono aachrono;
    Chrono aahyperchrono;
    Chrono totchrono;

    int profilecompmovenrep;
    int profilemulmovenrep;
    int simplehypermovenrep;
    int aishypermovenrep;

    double concaistuning;
    int concaisnstep;

    double acca1,acca2,acca3,acca4;
    double tota1,tota2,tota3,tota4;
    double accb1,accb2,accb3,accb4;
    double totb1,totb2,totb3,totb4;
    double centeracca1,centeracca2,centeracca3;
    double centertota1,centertota2,centertota3;
    double concacca1,concacca2,concacca3;
    double conctota1,conctota2,conctota3;
    double centeraccb1,centeraccb2,centeraccb3;
    double centertotb1,centertotb2,centertotb3;
    double concaccb1,concaccb2,concaccb3;
    double conctotb1,conctotb2,conctotb3;

	public:

    //-------------------
    // Construction and allocation
    // ------------------

	AAMutSelOmegaModel(string datafile, string treefile, int inaisnrep, double inaistuning, int inaisnstep) : aahypersuffstat(20) {

		data = new FileSequenceAlignment(datafile);
		codondata = new CodonSequenceAlignment(data, true);

		Nsite = codondata->GetNsite();    // # columns
		Ntaxa = codondata->GetNtaxa();

		std::cerr << "-- Number of sites: " << Nsite << std::endl;

		taxonset = codondata->GetTaxonSet();

		// get tree from file (newick format)
		tree = new Tree(treefile);

		// check whether tree and data fits together
		tree->RegisterWith(taxonset);

		tree->SetIndices();
		Nbranch = tree->GetNbranch();

		// Allocate();
        totchrono.Reset();
        aachrono.Reset();
        aahyperchrono.Reset();

        profilecompmovenrep = 5;
        profilemulmovenrep = 5;
        simplehypermovenrep = 100;
        aishypermovenrep = inaisnrep;
        concaistuning = inaistuning;
        concaisnstep = inaisnstep;

        acca1=acca2=acca3=acca4=0;
        tota1=tota2=tota3=tota4=0;
        accb1=accb2=accb3=accb4=0;
        totb1=totb2=totb3=totb4=0;

        centeracca1=centeracca2=centeracca3=0;
        centertota1=centertota2=centertota3=0;
        concacca1=concacca2=concacca3=0;
        conctota1=conctota2=conctota3=0;

        centeraccb1=centeraccb2=centeraccb3=0;
        centertotb1=centertotb2=centertotb3=0;
        concaccb1=concaccb2=concaccb3=0;
        conctotb1=conctotb2=conctotb3=0;
	}

    void Unfold()   {

		cerr << "-- unfold\n";
		phyloprocess->Unfold();
		cerr << phyloprocess->GetLogLikelihood() << '\n';
		std::cerr << "-- mapping substitutions\n";
		phyloprocess->ResampleSub();
		// Trace(cerr);
    }

	void Allocate()	{

		lambda = 10;
		branchlength = new BranchIIDGamma(*tree,1.0,lambda);
		lengthpathsuffstatarray = new PoissonSuffStatBranchArray(*tree);

		nucrelrate.assign(Nrr,0);
        Random::DirichletSample(nucrelrate,vector<double>(Nrr,1.0/Nrr),((double) Nrr));

		nucstat.assign(Nnuc,0);
        Random::DirichletSample(nucstat,vector<double>(Nnuc,1.0/Nnuc),((double) Nnuc));

		nucmatrix = new GTRSubMatrix(Nnuc,nucrelrate,nucstat,true);

        aacenter.assign(Naa,1.0/Naa);
        aaconc = ((double) Naa);
        aafitnessarray = new IIDDirichlet(Nsite,aacenter,aaconc);
        bkaafitnessarray = new IIDDirichlet(Nsite,aacenter,aaconc);

        omegahypermean = 1.0;
        omegahyperinvshape = 1.0;
		omega = 1.0;

        codonmatrixarray = new AAMutSelOmegaCodonSubMatrixArray(GetCodonStateSpace(), nucmatrix, aafitnessarray, omega);
		phyloprocess = new PhyloProcess(tree,codondata,branchlength,0,codonmatrixarray);
		pathsuffstatarray = new PathSuffStatArray(Nsite);
	}

    //-------------------
    // Accessors
    // ------------------

	CodonStateSpace* GetCodonStateSpace() const {
		return (CodonStateSpace*) codondata->GetStateSpace();
	}

    double GetOmega() const {
        return omega;
    }

    const PoissonSuffStatBranchArray* GetLengthPathSuffStatArray() const {
        return lengthpathsuffstatarray;
    }

    //-------------------
    // Setting and updating
    // ------------------

    void SetBranchLengths(const BranchSelector<double>& inbranchlength)    {
        branchlength->Copy(inbranchlength);
    }

    void SetOmega(double inomega)   {
        omega = inomega;
        UpdateCodonMatrices();
    }

    void SetOmegaHyperParameters(double inomegahypermean, double inomegahyperinvshape)   {
        omegahypermean = inomegahypermean;
        omegahyperinvshape = inomegahyperinvshape;
    }

    void SetNucRates(const std::vector<double>& innucrelrate, const std::vector<double>& innucstat) {
        nucrelrate = innucrelrate;
        nucstat = innucstat;
        UpdateMatrices();
    }

	void UpdateNucMatrix()	{
		nucmatrix->CopyStationary(nucstat);
		nucmatrix->CorruptMatrix();
	}

	void UpdateCodonMatrices()	{
        codonmatrixarray->SetOmega(omega);
		codonmatrixarray->UpdateCodonMatrices();
	}

    void UpdateCodonMatrix(int site)    {
        (*codonmatrixarray)[site].CorruptMatrix();
    }
		
    void UpdateMatrices()   {
        UpdateNucMatrix();
        UpdateCodonMatrices();
    }

    void NoUpdate() {}

    //-------------------
    // Priors and likelihood
    //-------------------

    double GetLogPrior() const {
        double total = 0;
        total += BranchLengthsHyperLogPrior();
        total += BranchLengthsLogPrior();
        total += NucRatesLogPrior();
        total += AAHyperLogPrior();
        total += AALogPrior();
        total += OmegaLogPrior();
        return total;
    }

	double GetLogLikelihood() const {
		return phyloprocess->GetLogLikelihood();
	}

    double GetLogProb() const   {
        return GetLogPrior() + GetLogLikelihood();
    }

	double BranchLengthsHyperLogPrior() const {
		return -lambda / 10;
	}

	double BranchLengthsLogPrior() const {
		return branchlength->GetLogProb();
	}

	// exponential of mean 1
	double OmegaLogPrior() const {
        double alpha = 1.0 / omegahyperinvshape;
        double beta = alpha / omegahypermean;
		return alpha * log(beta) - Random::logGamma(alpha) + (alpha-1) * log(omega) - beta*omega;
	}

    double AAHyperLogPrior() const {
        return Random::logGammaDensity(aaconc,20.0,1.0);
    }

    double AALogPrior() const {
        return aafitnessarray->GetLogProb();
    }

    double AALogPrior(int i) const {
        return aafitnessarray->GetLogProb(i);
    }

    double NucRatesLogPrior() const {
        return 0;
    }

    //-------------------
    // Suff Stat and suffstatlogprobs
    //-------------------

	double PathSuffStatLogProb() const {
		return pathsuffstatarray->GetLogProb(*codonmatrixarray);
	}

    double PathSuffStatLogProb(int site) const {
        return pathsuffstatarray->GetVal(site).GetLogProb(codonmatrixarray->GetVal(site));
    }

	double BranchLengthsHyperSuffStatLogProb() const {
		return lambdasuffstat.GetLogProb(1.0,lambda);
	}

    double AAHyperSuffStatLogProb() const   {
        return aahypersuffstat.GetLogProb(aacenter,aaconc);
    }

    //-------------------
    //  Log probs for MH moves
    //-------------------

    // for moving branch lengths hyperparameter lambda
    double BranchLengthsHyperLogProb() const {
        return BranchLengthsHyperLogPrior() + BranchLengthsHyperSuffStatLogProb();
    }

    // for moving nuc rates
    double NucRatesLogProb() const {
        return NucRatesLogPrior() + PathSuffStatLogProb();
    }

    // for moving aa hyper params (aacenter and aaconc)
    double AAHyperLogProb() const   {
        return AAHyperLogPrior() + AAHyperSuffStatLogProb();
    }

    //-------------------
    //  Moves 
    //-------------------

	double Move()	{
        ResampleSub(1.0);
        MoveParameters(30);
        return 1.0;
    }

    void ResampleSub(double frac)   {
        UpdateMatrices();
		phyloprocess->Move(frac);
    }

    void MoveParameters(int nrep)   {
		for (int rep=0; rep<nrep; rep++)	{

            totchrono.Start();
			ResampleBranchLengths();
			MoveBranchLengthsHyperParameter();

			CollectPathSuffStat();

			MoveNucRates();
			MoveOmega();

            aachrono.Start();
            CompMoveAA(profilecompmovenrep);
            MulMoveAA(profilemulmovenrep);
            aachrono.Stop();

            aahyperchrono.Start();
            MoveAAHyperParameters(simplehypermovenrep);
            AISMoveAAHyperParameters(aishypermovenrep);
            aahyperchrono.Stop();

            totchrono.Stop();
		}
	}

    void CollectLengthSuffStat()    {
		lengthpathsuffstatarray->Clear();
        lengthpathsuffstatarray->AddLengthPathSuffStat(*phyloprocess);
    }

	void ResampleBranchLengths()	{
        CollectLengthSuffStat();
		branchlength->GibbsResample(*lengthpathsuffstatarray);
	}

	void MoveBranchLengthsHyperParameter()	{

		lambdasuffstat.Clear();
		branchlength->AddSuffStat(lambdasuffstat);
        ScalingMove(lambda,1.0,10,&AAMutSelOmegaModel::BranchLengthsHyperLogProb,&AAMutSelOmegaModel::NoUpdate,this);
        ScalingMove(lambda,0.3,10,&AAMutSelOmegaModel::BranchLengthsHyperLogProb,&AAMutSelOmegaModel::NoUpdate,this);
		branchlength->SetScale(lambda);
	}

	void CollectPathSuffStat()	{

		pathsuffstatarray->Clear();
        pathsuffstatarray->AddSuffStat(*phyloprocess);
	}

	void MoveOmega()	{

		omegasuffstat.Clear();
		omegasuffstat.AddSuffStat(*codonmatrixarray,*pathsuffstatarray);
        double alpha = 1.0 / omegahyperinvshape;
        double beta = alpha / omegahypermean;
		omega = Random::GammaSample(alpha + omegasuffstat.GetCount(), beta + omegasuffstat.GetBeta());
		UpdateCodonMatrices();
	}

	void MoveNucRates()	{

        // UpdateMatrices();

        ProfileMove(nucrelrate,0.1,1,3,&AAMutSelOmegaModel::NucRatesLogProb,&AAMutSelOmegaModel::UpdateMatrices,this);
        ProfileMove(nucrelrate,0.03,3,3,&AAMutSelOmegaModel::NucRatesLogProb,&AAMutSelOmegaModel::UpdateMatrices,this);
        ProfileMove(nucrelrate,0.01,3,3,&AAMutSelOmegaModel::NucRatesLogProb,&AAMutSelOmegaModel::UpdateMatrices,this);

        ProfileMove(nucstat,0.1,1,3,&AAMutSelOmegaModel::NucRatesLogProb,&AAMutSelOmegaModel::UpdateMatrices,this);
        ProfileMove(nucstat,0.01,1,3,&AAMutSelOmegaModel::NucRatesLogProb,&AAMutSelOmegaModel::UpdateMatrices,this);

        // UpdateMatrices();
	}

    double AISMoveAAHyperConcentration(double tuning, int nstep)  {

        // backup all fitness profiles
        bkaafitnessarray->Copy(*aafitnessarray);

        double aaconc1 = aaconc;
        double m = tuning*(Random::Uniform() - 0.5);
        double e = exp(m);
        double aaconc2 = aaconc * e;
        double loghastings = m;

        double logratio = 0;

        for (int i=0; i<nstep; i++) {

            aahypersuffstat.Clear();
            aafitnessarray->AddSuffStat(aahypersuffstat);

            aaconc = ((nstep-i)*aaconc1 + i*aaconc2)/nstep;
            logratio -= AAHyperSuffStatLogProb();

            aaconc = ((nstep-i-1)*aaconc1 + (i+1)*aaconc2)/nstep;
            logratio += AAHyperSuffStatLogProb();

            if (i < nstep-1)    {
                aafitnessarray->SetConcentration(aaconc);
                if (Random::Uniform() < 0.5)    {
                    MoveAAGamma(3.0,1);
                }
                else    {
                    MoveAAGamma(1.0,1);
                }
                /*
                MulMoveAA(1);
                CompMoveAA(1);
                */
            }
        }

        aaconc = aaconc1;
        logratio -= AAHyperLogPrior();

        aaconc = aaconc2;
        logratio += AAHyperLogPrior();

        logratio += loghastings;

        int accepted = (log(Random::Uniform()) < logratio);

        if (! accepted) {
            aaconc = aaconc1;
            // restore fitness profiles
            aafitnessarray->Copy(*bkaafitnessarray);
        }
        else    {
            aaconc = aaconc2;
        }
        aafitnessarray->SetConcentration(aaconc);
        return accepted;
    }

    double AISMoveAAHyperCenter(double tuning, int n, int nstep)  {

        // backup all fitness profiles
        bkaafitnessarray->Copy(*aafitnessarray);

        // propose new value
        vector<double> aacenter1 = aacenter;
        double loghastings = Random::ProfileProposeMove(aacenter,Naa,tuning,n);
        vector<double> aacenter2 = aacenter;

        double logratio = 0;

        for (int i=0; i<nstep; i++) {

            aahypersuffstat.Clear();
            aafitnessarray->AddSuffStat(aahypersuffstat);

            for (int k=0; k<Naa; k++)   {
                aacenter[k] = ((nstep-i)*aacenter1[k] + i*aacenter2[k])/nstep;
            }
            logratio -= AAHyperSuffStatLogProb();

            for (int k=0; k<Naa; k++)   {
                aacenter[k] = ((nstep-i-1)*aacenter1[k] + (i+1)*aacenter2[k])/nstep;
            }
            logratio += AAHyperSuffStatLogProb();

            if (i < nstep-1)    {
                aafitnessarray->SetCenter(aacenter);
                MulMoveAA(1);
                CompMoveAA(1);
            }
        }

        aacenter = aacenter1;
        logratio -= AAHyperLogPrior();

        aacenter = aacenter2;
        logratio += AAHyperLogPrior();

        logratio += loghastings;

        int accepted = (log(Random::Uniform()) < logratio);

        if (! accepted) {
            aacenter = aacenter1;
            // restore fitness profiles
            aafitnessarray->Copy(*bkaafitnessarray);
        }
        else    {
            aacenter = aacenter2;
        }
        aafitnessarray->SetCenter(aacenter);
        return accepted;
    }

    void AISMoveAAHyperParameters(int nrep)   {

        for (int rep=0; rep<nrep; rep++)  {
            // AISMoveAAHyperCenter(0.1,3,10);
            concaccb1 += AISMoveAAHyperConcentration(concaistuning,concaisnstep);
            conctotb1++;
        }
    }

    void MoveAAHyperParameters(int nrep)    {
        aahypersuffstat.Clear();
        aafitnessarray->AddSuffStat(aahypersuffstat);
        for (int rep=0; rep<nrep; rep++)  {
            centeracca1 += ProfileMove(aacenter,0.1,1,10,&AAMutSelOmegaModel::AAHyperLogProb,&AAMutSelOmegaModel::NoUpdate,this);
            centertota1++;
            centeracca2 += ProfileMove(aacenter,0.03,3,10,&AAMutSelOmegaModel::AAHyperLogProb,&AAMutSelOmegaModel::NoUpdate,this);
            centertota2++;
            centeracca3 += ProfileMove(aacenter,0.01,3,10,&AAMutSelOmegaModel::AAHyperLogProb,&AAMutSelOmegaModel::NoUpdate,this);
            centertota3++;
            concacca1 += ScalingMove(aaconc,1.0,10,&AAMutSelOmegaModel::AAHyperLogProb,&AAMutSelOmegaModel::NoUpdate,this);
            conctota1++;
            concacca2 += ScalingMove(aaconc,0.3,10,&AAMutSelOmegaModel::AAHyperLogProb,&AAMutSelOmegaModel::NoUpdate,this);
            conctota2++;
            concacca3 += ScalingMove(aaconc,0.1,10,&AAMutSelOmegaModel::AAHyperLogProb,&AAMutSelOmegaModel::NoUpdate,this);
            conctota3++;
        }
        aafitnessarray->SetCenter(aacenter);
        aafitnessarray->SetConcentration(aaconc);
    }

    double CompMoveAA(int nrep) {
        accb1 += MoveAA(1.0,1,nrep);
        accb2 += MoveAA(1.0,3,nrep);
        accb3 += MoveAA(0.3,3,nrep);
        accb4 += MoveAA(0.1,3,nrep);
        totb1++;
        totb2++;
        totb3++;
        totb4++;
        return 1.0;
    }

    double MulMoveAA(int nrep) {
        acca1 += MoveAAGamma(3.0,nrep);
        acca2 += MoveAAGamma(1.0,nrep);
        acca3 += MoveAAGamma(0.3,nrep);
        acca4 += MoveAAGamma(0.1,nrep);
        tota1++;
        tota2++;
        tota3++;
        tota4++;
        return 1.0;
    }

	double MoveAA(double tuning, int n, int nrep)	{
		double nacc = 0;
		double ntot = 0;
        vector<double> bk(20,0);
        for (int i=0; i<Nsite; i++) {
            vector<double>& aa = (*aafitnessarray)[i];
            for (int rep=0; rep<nrep; rep++)	{
                bk = aa;
                double deltalogprob = -AALogPrior(i) - PathSuffStatLogProb(i);
                double loghastings = Random::ProfileProposeMove(aa,Naa,tuning,n);
                deltalogprob += loghastings;
                UpdateCodonMatrix(i);
                deltalogprob += AALogPrior(i) + PathSuffStatLogProb(i);
                int accepted = (log(Random::Uniform()) < deltalogprob);
                if (accepted)	{
                    nacc ++;
                }
                else	{
                    aa = bk;
                    UpdateCodonMatrix(i);
                }
                ntot++;
            }
        }
		return nacc/ntot;
	}

    double GammaAALogPrior(const vector<double>& x) {
        double total = 0;
        for (int l=0; l<Naa; l++)   {
            total += (aaconc*aacenter[l] -1)*log(x[l]) - x[l] - Random::logGamma(aaconc*aacenter[l]);
        }
        return total;
    }

	double MoveAAGamma(double tuning, int nrep)	{

		double nacc = 0;
		double ntot = 0;
        for (int i=0; i<Nsite; i++) {

            vector<double>& aa = (*aafitnessarray)[i];
            vector<double> x(Naa,0);
            double z = Random::sGamma(aaconc);
            for (int l=0; l<Naa; l++)   {
                x[l] = z*aa[l];
            }

            double bkz = z;
            vector<double> bkx = x;
            vector<double> bkaa = aa;

            for (int rep=0; rep<nrep; rep++)	{

                double deltalogprob = -GammaAALogPrior(x) - PathSuffStatLogProb(i);

                double loghastings = 0;
                z = 0;
                for (int l=0; l<Naa; l++)   {
                    double m = tuning * (Random::Uniform() - 0.5);
                    double e = exp(m);
                    x[l] *= e;
                    z += x[l];
                    loghastings += m;
                }
                for (int l=0; l<Naa; l++)   {
                    aa[l] = x[l]/z;
                }

                deltalogprob += loghastings;

                UpdateCodonMatrix(i);

                deltalogprob += GammaAALogPrior(x) + PathSuffStatLogProb(i);
                int accepted = (log(Random::Uniform()) < deltalogprob);
                if (accepted)	{
                    nacc ++;
                    bkaa = aa;
                    bkx = x;
                    bkz = z;
                }
                else	{
                    aa = bkaa;
                    x = bkx;
                    z = bkz;
                    UpdateCodonMatrix(i);
                }
                ntot++;
            }
        }
		return nacc/ntot;
	}

    //-------------------
    // Traces and Monitors
    // ------------------

	void TraceHeader(std::ostream& os) const {
		os << "#logprior\tlnL\tlength\t";
		os << "omega\t";
        os << "aaent\t";
        os << "aaconc\t";
        os << "aacenter\t";
		os << "statent\t";
		os << "rrent\n";
	}

	void Trace(ostream& os) const {	
		os << GetLogPrior() << '\t';
		os << GetLogLikelihood() << '\t';
        // 3x: per coding site (and not per nucleotide site)
        os << 3*branchlength->GetTotalLength() << '\t';
		os << omega << '\t';
        os << aafitnessarray->GetMeanEntropy() << '\t';
        os << aaconc << '\t';
        os << Random::GetEntropy(aacenter) << '\t';
		os << Random::GetEntropy(nucstat) << '\t';
		os << Random::GetEntropy(nucrelrate) << '\n';

    }

	void Monitor(ostream& os) const {
        os << "prop time in aa moves: " << aachrono.GetTime() / totchrono.GetTime() << '\n';
        os << "prop time in aa hyper moves: " << aahyperchrono.GetTime() / totchrono.GetTime() << '\n';

        if (profilemulmovenrep) {
            os << acca1/tota1 << '\n';
            os << acca2/tota2 << '\n';
            os << acca3/tota3 << '\n';
            os << acca4/tota4 << '\n';
        }
        os << '\n';

        if (profilecompmovenrep)  {
            os << accb1/totb1 << '\n';
            os << accb2/totb2 << '\n';
            os << accb3/totb3 << '\n';
            os << accb4/totb4 << '\n';
        }
        os << '\n';

        if (simplehypermovenrep)    {
            os << centeracca1/centertota1 << '\n';
            os << centeracca2/centertota2 << '\n';
            os << centeracca3/centertota3 << '\n';
            os << concacca1/conctota1 << '\n';
            os << concacca2/conctota2 << '\n';
            os << concacca3/conctota3 << '\n';
        }
        os << '\n';

        if (aishypermovenrep)   {
            os << "ais: " << aishypermovenrep << '\t' << concaistuning << '\t' << concaisnstep << '\n';
            os << concaccb1/conctotb1 << '\n';
        }
        os << '\n';
    }

	void FromStream(istream& is) {}
	void ToStream(ostream& os) const {}

};


