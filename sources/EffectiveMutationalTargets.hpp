

// take 
//  a multiple sequence alignment, 
//  a tree,
//  an array of branch codon matrices (BranchSelector)
// for each leaf: 
//  compute empirical freqs in corresponding species in alignment
//  compute number of targets
//  print out

void RecursiveAddEffectiveMutationalTargets(const Link* from, const CodonSequenceAlignment& ali, const BranchSelector<MGOmegaCodonSubMatrix>& mats, std::vector<double>& syn, std::vector<double>& nonsyn)    {

    if (from->isLeaf()) {
        int tax = from->GetNode()->GetIndex();
        int branch = from->GetBranch()->GetIndex();
        const MGOmegaCodonSubMatrix& mat = mats.GetVal(branch);
        std::vector<int> counts = ali.GetEmpiricalCounts(tax);
        mat.GetMeanEffectiveMutationalTargets(counts, syn[tax], nonsyn[tax]);
    }
    for (const Link* link=from->Next(); link!=from; link=link->Next())  {
        RecursiveAddEffectiveMutationalTargets(link->Out(), ali, mats, syn, nonsyn);
    }
}


