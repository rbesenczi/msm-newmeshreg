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

        const int status_bar_width = 70;

        for(int i = 0; i < mc_iterations; ++i)
        {
            for(int node = 0; node < energy->getNumNodes(); ++node)
            {
                double old_unary_cost = 0.0, new_unary_cost = 0.0,
                       old_triplet_cost = 0.0, new_triplet_cost = 0.0;

                int label = distribution(gen);

                old_unary_cost = energy->computeUnaryCost(node, labeling[node]);
                new_unary_cost = energy->computeUnaryCost(node, label);

                for(int ntr = 0; ntr < cp_grid.get_total_triangles(node); ++ntr)
                {
                    const newresampler::Triangle& tr = cp_grid.get_triangle_from_vertex(node, ntr);
                    old_triplet_cost += energy->computeTripletCostTri(tr.get_no(), labeling[tr.get_vertex_no(0)], labeling[tr.get_vertex_no(1)], labeling[tr.get_vertex_no(2)]);
                    new_triplet_cost += energy->computeTripletCostTri(tr.get_no(), label, labeling[tr.get_vertex_no(1)], labeling[tr.get_vertex_no(2)]);
                }

                if (new_unary_cost + new_triplet_cost < old_unary_cost + old_triplet_cost)
                    labeling[node] = label;
            }
            //--
            // The following few lines are for a status bar, don't bother with them...
            std::cout << "Monte Carlo optimisation progress [";
            double progress = (double)i/mc_iterations;
            int pos = status_bar_width * progress;
            for (int k = 0; k < status_bar_width; ++k)
            {
                if (k < pos) std::cout << "=";
                else if (k == pos) std::cout << ">";
                else std::cout << " ";
            }
            std::cout << "] " << std::ceil(progress * 100.0) << " %\r";
            std::cout.flush();
            //--
        }
        std::cout << std::endl;
        return energy->evaluateTotalCostSum();
    }
};

} //namespace newmeshreg

#endif //NEWMESHREG_MCMC_OPT_H
