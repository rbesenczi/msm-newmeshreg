#ifndef NEWMESHREG_MCMC_OPT_H
#define NEWMESHREG_MCMC_OPT_H

#include "DiscreteModel.h"

namespace newmeshreg {

class MCMC {
public:
    static double optimise(const std::shared_ptr<DiscreteModel>& energy, bool verbose, int numthreads) {

        const int *triplets = energy->getTriplets();
        const int *pairs = energy->getPairs();
        int *labeling = energy->getLabeling();
        auto costfnct = energy->getCostFunction();
        double* unary_costs = costfnct->getUnaryCosts();
        double* pairwise_costs = costfnct->getPairwiseCosts();
        double* triplet_costs = costfnct->getTripletCosts();
        const int mc_iterations = 1000;

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> distribution(0, energy->getNumLabels()-1);

        for(int i = 0; i < mc_iterations; ++i) {
            for(int node = 0; node < energy->getNumNodes(); ++node) {
                double current_energy = unary_costs[node+labeling[node]] + pairwise_costs[node*2+labeling[node]] + triplet_costs[node*3+labeling[node]];
                int label = distribution(gen);
                double new_energy = unary_costs[node+label] + pairwise_costs[node*2+label] + triplet_costs[node*3+label];
                if(std::abs(new_energy) < std::abs(current_energy))
                {
                    std::cout << "\tLabel changed..." << "\tnode==" << node << "\toldlabel==" << labeling[node] << "\tnewlabel==" << label << "\tcurrent_energy==" << current_energy << "\tnew_energy==" << new_energy << std::endl;
                    labeling[node] = label;
                }
            }
        }

        return energy->evaluateTotalCostSum();
    }
};

} //namespace newmeshreg

#endif //NEWMESHREG_MCMC_OPT_H
