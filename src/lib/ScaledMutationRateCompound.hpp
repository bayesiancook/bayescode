#pragma once

#include "MultivariateProcess.hpp"
#include "ScaledMutationRate.hpp"

/**
 * \brief A ScaledMutationRate that returns the scaled mutation rate (theta=4*Ne*u) as a compound
 * parameter given Ne (population size) and u (mutation rate per generation) of a given taxon.
 *
 * Useful for implementing models assuming each taxon has its own Ne (population size) and u
 * (mutation rate per generation). Moreover Ne and u are assumed to follow a brownian process
 * (NodeProcess).
 */
class NodeProcessScaledMutationRate : public ScaledMutationRate {
  public:
    //! \brief Constructor, taking as arguments the Nbr of taxa, the scaling factor, the mutation
    //! rate (a NodeProcess) and the population size (a NodeProcess).
    NodeProcessScaledMutationRate(double const &intheta_scale,
        NodeProcess const *innode_rates, NodeProcess const *innode_popsize, TaxonSet const &taxon);

    ~NodeProcessScaledMutationRate() override = default;

    double GetTheta(int taxon) const override;

    double &operator[](int taxon) { return theta[taxon]; }

    void Update();

  private:
    double const &theta_scale;
    NodeProcess const *node_rates;
    NodeProcess const *node_popsize;
    std::vector<double> theta;

    std::vector<Tree::NodeIndex> reverse_index_table;
};

class BranchWiseProcessScaledMutationRate : public ScaledMutationRate {
  public:
    //! \brief Constructor, taking as arguments the Nbr of taxa, the scaling factor, the mutation
    //! rate (a NodeProcess) and the population size (a NodeProcess).
    BranchWiseProcessScaledMutationRate(double const &intheta_scale,
        LeafMultivariateProcess &inleaf_multivariate_process, TaxonSet const &taxon);

    ~BranchWiseProcessScaledMutationRate() override = default;

    double GetTheta(int taxon) const override { return theta[taxon]; }

    double &operator[](int taxon) { return theta[taxon]; }

    void SlidingTaxonMove(int taxon, int dimension, double m);

    void Update();

  private:
    double const &theta_scale;
    LeafMultivariateProcess &leaf_multivariate_process;
    std::vector<double> theta;

    std::vector<Tree::NodeIndex> reverse_index_table;
};