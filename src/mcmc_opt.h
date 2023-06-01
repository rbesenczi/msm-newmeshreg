#ifndef NEWMESHREG_MCMC_OPT_H
#define NEWMESHREG_MCMC_OPT_H

#include "DiscreteModel.h"

namespace newmeshreg {

class MCMC {
public:
    static double optimise(const std::shared_ptr<NonLinearSRegDiscreteModel>& energy, bool verbose, int numthreads) {

        const int mc_iterations = 100;

        int* labeling = energy->getLabeling();
        const auto cp_grid = energy->get_CPgrid();

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> distribution(0, energy->getNumLabels()-1);

        for(int i = 0; i < mc_iterations; ++i)
        {
            for(int node = 0; node < energy->getNumNodes(); ++node)
            {
                int label = distribution(gen);
                double old_unary_cost = energy->computeUnaryCost(node, labeling[node]);
                double new_unary_cost = energy->computeUnaryCost(node, label);

                double old_pairwise_cost = 0.0;
                double new_pairwise_cost = 0.0;

                double old_triplet_cost = 0.0;
                double new_triplet_cost = 0.0;

                int num_neighbours = 0;

                for(auto it = cp_grid.nbegin(node); it != cp_grid.nend(node); it++)
                {
                    old_pairwise_cost += energy->computePairwiseCost(node + num_neighbours*2,
                                                                    labeling[node],
                                                                    labeling[*it]);
                    new_pairwise_cost += energy->computePairwiseCost(node + num_neighbours*2,
                                                                    label,
                                                                    labeling[*it]);
                    num_neighbours++;
                }

                //compute triplet cost.

                if(old_unary_cost+old_pairwise_cost+old_triplet_cost > new_unary_cost+new_pairwise_cost+new_triplet_cost)
                {
                    labeling[node] = label;
                    if(verbose)
                        std::cout << "Node " << node << " assigned to label " << label << ".\n"
                                  << "\told_unary_cost==" << old_unary_cost << "\tnew_unary_cost==" << new_unary_cost <<
                                  "\told_pairwise_cost==" << old_pairwise_cost << "\tnew_pairwise_cost==" << new_pairwise_cost << std::endl;
                }
            }
        }
        return energy->evaluateTotalCostSum();
    }
};

} //namespace newmeshreg

#endif //NEWMESHREG_MCMC_OPT_H
