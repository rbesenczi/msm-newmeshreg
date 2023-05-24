#ifndef NEWMESHREG_MCMC_OPT_H
#define NEWMESHREG_MCMC_OPT_H

#include "DiscreteModel.h"

namespace newmeshreg {

class MCMC {
    static double optimise(const std::shared_ptr<DiscreteModel>& energy, bool verbose, int numthreads) {

        const int *triplets = energy->getTriplets();
        const int *pairs = energy->getPairs();
        int *labeling = energy->getLabeling();
        const int mc_iteration = 1000;

        double initEnergy = energy->evaluateTotalCostSum();
        double lastEnergy = initEnergy;

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> distribution(0, energy->getNumLabels()-1);

        for(int i = 0; i < mc_iteration; ++i) {
            int label = distribution(gen);
        }

        lastEnergy = energy->evaluateTotalCostSum();

        return lastEnergy;
    }
};

} //namespace newmeshreg

#endif //NEWMESHREG_MCMC_OPT_H
