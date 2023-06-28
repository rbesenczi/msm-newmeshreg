#ifndef NEWMESHREG_MCMC_OPT_H
#define NEWMESHREG_MCMC_OPT_H

#include "DiscreteModel.h"
#include <algorithm>

namespace newmeshreg {

class MCMC {
public:
    static double optimise(const std::shared_ptr<NonLinearSRegDiscreteModel>& energy, bool verbose, int mciters) {

        int* labeling = energy->getLabeling();
        const int num_nodes = energy->getNumNodes();
        const double* unary_costs = energy->getCostFunction()->getUnaryCosts();
        const int* triplets = energy->getTriplets();

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> distribution(0, energy->getNumLabels()-1);

        for(int i = 0; i < mciters; ++i)
        {
            for(int triplet = 0; triplet < energy->getNumTriplets(); ++triplet)
            {
                std::vector<double> triplet_data(4);

                int label = distribution(gen);

                const int nodeA = triplets[triplet*3];
                const int nodeB = triplets[triplet*3+1];
                const int nodeC = triplets[triplet*3+2];

                triplet_data[0] = energy->computeTripletCost(triplet,labeling[nodeA],labeling[nodeB],labeling[nodeC]) +
                                  unary_costs[labeling[nodeC] * num_nodes + nodeC] + unary_costs[labeling[nodeB] * num_nodes + nodeB] + unary_costs[labeling[nodeA] * num_nodes + nodeA];
                triplet_data[1] = energy->computeTripletCost(triplet,labeling[nodeA],labeling[nodeB],label) +
                                  unary_costs[label * num_nodes + nodeC] + unary_costs[labeling[nodeB] * num_nodes + nodeB] + unary_costs[labeling[nodeA] * num_nodes + nodeA];
                triplet_data[2] = energy->computeTripletCost(triplet,labeling[nodeA],label,labeling[nodeC]) +
                                  unary_costs[label * num_nodes + nodeB] + unary_costs[labeling[nodeC] * num_nodes + nodeC] + unary_costs[labeling[nodeA] * num_nodes + nodeA];
                triplet_data[3] = energy->computeTripletCost(triplet,label,labeling[nodeB],labeling[nodeC]) +
                                  unary_costs[label * num_nodes + nodeA] + unary_costs[labeling[nodeC] * num_nodes + nodeC] + unary_costs[labeling[nodeB] * num_nodes + nodeB];

                int index = (int)std::ranges::distance(triplet_data.begin(), std::ranges::min_element(triplet_data));

                if(index != 0) {
                    switch (index) {
                        case 1:
                            labeling[nodeC] = label;
                            break;
                        case 2:
                            labeling[nodeB] = label;
                            break;
                        case 3:
                            labeling[nodeA] = label;
                            break;
                    }
                }
            }

            if (verbose && i % 100 == 0)
            {
                std::cout << "MC iter " << i << '/' << mciters << '\t';
                energy->evaluateTotalCostSum();
            }
        }

        if (verbose)
            std::cout << "MC iter " << mciters << '/' << mciters << '\t';

        return energy->evaluateTotalCostSum();
    }
};

} //namespace newmeshreg

#endif //NEWMESHREG_MCMC_OPT_H
