#ifndef NEWMESHREG_MCMC_OPT_H
#define NEWMESHREG_MCMC_OPT_H

#include "DiscreteModel.h"

namespace newmeshreg {

class MCMC {
public:
    static double optimise(const std::shared_ptr<NonLinearSRegDiscreteModel>& energy, bool verbose, int mciters) {

        int* labeling = energy->getLabeling();
        const int* triplets = energy->getTriplets();

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> distribution(0, energy->getNumLabels()-1);

        const int status_bar_width = 70;

        energy->evaluateTotalCostSum();

        for(int i = 0; i < mciters; ++i)
        {
            for(int triplet = 0; triplet < energy->getNumTriplets(); ++triplet)
            {
                std::vector<double> triplet_data(4);

                int label = distribution(gen);

                const int nodeA = triplets[triplet*3];
                const int nodeB = triplets[triplet*3+1];
                const int nodeC = triplets[triplet*3+2];

                triplet_data[0] = energy->computeTripletCost(triplet,labeling[nodeA],labeling[nodeB],labeling[nodeC]);	//000
                triplet_data[1] = energy->computeTripletCost(triplet,labeling[nodeA],labeling[nodeB],label);			//001
                triplet_data[2] = energy->computeTripletCost(triplet,labeling[nodeA],label,labeling[nodeC]);			//010
                triplet_data[3] = energy->computeTripletCost(triplet,label,labeling[nodeB],labeling[nodeC]);			//100

                auto minimum = std::min_element(triplet_data.begin(), triplet_data.end());
                int num = 0;
                for(auto it = triplet_data.begin(); it != minimum; it++)
                    num++;

                if(num > 0) {
                    switch (num) {
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

            if(verbose) {
                // The following few lines are for a status bar, don't bother with them...
                std::cout << "MC progress [";
                double progress = (double) i / mciters;
                int pos = status_bar_width * progress;
                for (int k = 0; k < status_bar_width; ++k) {
                    if (k < pos) std::cout << "=";
                    else if (k == pos) std::cout << ">";
                    else std::cout << " ";
                }
                std::cout << "] " << std::ceil(progress * 100.0) << " %\r";
                std::cout.flush();
            }

        }
        std::cout << std::endl;

        return energy->evaluateTotalCostSum();
    }
};

} //namespace newmeshreg

#endif //NEWMESHREG_MCMC_OPT_H
