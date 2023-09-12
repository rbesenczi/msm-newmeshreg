#ifndef NEWMESHREG_MCMC_OPT_H
#define NEWMESHREG_MCMC_OPT_H

#include "DiscreteModel.h"

namespace newmeshreg {

class MCMC {
public:
    static double optimise(const std::shared_ptr<NonLinearSRegDiscreteModel>& energy, bool verbose, int mciters) {

        int* labeling = energy->getLabeling();
        const int num_nodes = energy->getNumNodes();
        const double* unary_costs = energy->getCostFunction()->getUnaryCosts();
        const int* triplets = energy->getTriplets();

        if(verbose) { std::cout << "Initial "; energy->evaluateTotalCostSum(); }

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> distribution(0, energy->getNumLabels()-1);

        for(int i = 0; i < mciters; ++i)
        {
            for(int triplet = 0; triplet < energy->getNumTriplets(); ++triplet)
            {
                int label = distribution(gen);

                const int nodeA = triplets[triplet*3];
                const int nodeB = triplets[triplet*3+1];
                const int nodeC = triplets[triplet*3+2];

                if ((energy->computeTripletCost(triplet, label, labeling[nodeB], labeling[nodeC]) +
                    unary_costs[label * num_nodes + nodeA]) <
                    (energy->computeTripletCost(triplet, labeling[nodeA], labeling[nodeB], labeling[nodeC]) +
                    unary_costs[labeling[nodeA] * num_nodes + nodeA]))
                        labeling[nodeA] = label;
            }

            if (verbose && i % 1000 == 0 && i > 0)
            {
                std::cout << "MC iter " << i << '/' << mciters << '\t';
                energy->evaluateTotalCostSum();
            }
        }

        if (verbose) std::cout << "MC iter " << mciters << '/' << mciters << '\t';

        return energy->evaluateTotalCostSum();
    }
};

} //namespace newmeshreg

#endif //NEWMESHREG_MCMC_OPT_H
