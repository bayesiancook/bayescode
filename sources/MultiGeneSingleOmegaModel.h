/*
simple model: gene-specific omegas, iid from a gamma distribution
shared nuc rates and branch lengths across all genes
each slave allocates a small number of gene-specific models
master deals with branch lengths, nuc stat and shape and scale params of gamma distribution of omega's across genes

move schedule:

starting:
master broadcasts global params / slaves receive global params and update phyloprocess

master collects branch length suff stat across genes / slaves compile BL suff stats and send them to master
master resamples branch lengths (and hyperparams) and broadcasts branch lengths / slaves receive branch lengths

master collects nuc suffstat across genes / slave compile nuc suff stats and send them to master
master resamples nuc rates and broadcasts them / slave receive and update

(1) master receive omegas, move alpha and beta and send alpha and beta / slave move omega | alpha, beta, then send omega, then receive new alpha/beta
(2) master collects omega suff stat across genes, move alpha beta (integrated), resample omega's and send them to slave

- use SingleOmegaModel as a single class, usable under all conditions (with flags indicating which aspects are dependent / independent)
*/

class MultiGeneSingleOmegaModel	{

	// variable indicating whether master (0) or slave (1--nproc)
	// simple MPI scheduler, with slave or master-specific behavior

};

