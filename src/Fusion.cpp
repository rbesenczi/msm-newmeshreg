#include "Fusion.h"

namespace newmeshreg {

template<typename OPTIMIZER>
void Fusion::reduce_and_convert(ELCReduce::PBF<REAL>& pbf, OPTIMIZER& MODEL, Reduction mode) {

    switch(mode)
    {
        case ELC_HOCR:
        {
            ELCReduce::PBF<REAL> qpbf;
            pbf.reduceHigher();
            pbf.toQuadratic(qpbf, pbf.maxID()+1);
            // Reduce the remaining higher-order terms using HOCR adding auxiliary variables
            qpbf.convert(MODEL, qpbf.maxID()+1);
            pbf.clear();
            qpbf.clear();
            break;
        }
        case ELC_APPROX:
        {
            pbf.reduceHigherApprox();
            pbf.convert(MODEL, pbf.maxID()+1);
            pbf.clear();
            break;
        }
        case HOCR:
        {
            ELCReduce::PBF<REAL> qpbf;
            pbf.toQuadratic(qpbf, pbf.maxID()+1); // Reduce to Quadratic pseudo-Boolean function using HOCR.
            qpbf.convert(MODEL, qpbf.maxID()+1);
            pbf.clear();
            qpbf.clear();
            break;
        }
    }
}

double Fusion::optimize(const std::shared_ptr<DiscreteModel>& energy, Reduction reductionMode, bool verbose, int numthreads) {

    const int *triplets = energy->getTriplets();
    int *labeling = energy->getLabeling();
    std::shared_ptr<DiscreteModelDummy> FPDMODEL = std::make_shared<DiscreteModelDummy>();

    double initEnergy = energy->evaluateTotalCostSum();
    double lastEnergy = initEnergy;
    double sumlabeldiff = 0.0;

    for(int sweep = 0; sweep < NUM_SWEEPS; ++sweep)
    {
        for(int label = 0; label < energy->getNumLabels(); ++label)
        {
            sumlabeldiff = 0;
            ELCReduce::PBF<REAL> pbf;
            int improveCounter = 0;
            double ratioUnlabeled = 0.0;
            int nodesChanged = 0;

            std::vector<UnaryData> unary_data(energy->getNumNodes());

            #pragma omp parallel for num_threads(numthreads)
            for(int node = 0; node < energy->getNumNodes(); ++node)
            {
                unary_data[node].buffer[0] = energy->computeUnaryCost(node,labeling[node]);	//0
                unary_data[node].buffer[1] = energy->computeUnaryCost(node,label);				//1
                #pragma omp critical
                sumlabeldiff += abs(label-labeling[node]);
            }

            if(sumlabeldiff > 0)
            {
                for (int node = 0; node < energy->getNumNodes(); ++node)
                    pbf.AddUnaryTerm(node, unary_data[node].buffer[0], unary_data[node].buffer[1]);

                std::vector<TripletData> triplet_data(energy->getNumTriplets());

                #pragma omp parallel for num_threads(numthreads)
                for (int triplet = 0; triplet < energy->getNumTriplets(); ++triplet)
                {
                    const int nodeA = triplets[triplet*3];
                    const int nodeB = triplets[triplet*3+1];
                    const int nodeC = triplets[triplet*3+2];

                    triplet_data[triplet].buffer[0] = energy->computeTripletCost(triplet,labeling[nodeA],labeling[nodeB],labeling[nodeC]);	//000
                    triplet_data[triplet].buffer[1] = energy->computeTripletCost(triplet,labeling[nodeA],labeling[nodeB],label);			//001
                    triplet_data[triplet].buffer[2] = energy->computeTripletCost(triplet,labeling[nodeA],label,labeling[nodeC]);			//010
                    triplet_data[triplet].buffer[3] = energy->computeTripletCost(triplet,labeling[nodeA],label,label);						//011
                    triplet_data[triplet].buffer[4] = energy->computeTripletCost(triplet,label,labeling[nodeB],labeling[nodeC]);			//100
                    triplet_data[triplet].buffer[5] = energy->computeTripletCost(triplet,label,labeling[nodeB],label);						//101
                    triplet_data[triplet].buffer[6] = energy->computeTripletCost(triplet,label,label,labeling[nodeC]);						//110
                    triplet_data[triplet].buffer[7] = energy->computeTripletCost(triplet,label,label,label);								//111
                }

                for (int triplet = 0; triplet < energy->getNumTriplets(); ++triplet)
                {
                    const int nodeA = triplets[triplet*3];
                    const int nodeB = triplets[triplet*3+1];
                    const int nodeC = triplets[triplet*3+2];
                    int node_ids[3] = { nodeA, nodeB, nodeC };
                    pbf.AddHigherTerm(3, node_ids, triplet_data[triplet].buffer);
                }

                FPDMODEL->reset();
                reduce_and_convert(pbf,*FPDMODEL, reductionMode);
                FPDMODEL->initialise();
                int* Labels = FPDMODEL->getLabeling();

                FPD::FastPD opt(FPDMODEL, 100 );
                double newEnergy = opt.run();
                opt.getLabeling(Labels);

                for(int node = 0; node < energy->getNumNodes(); ++node)
                {
                    if(labeling[node] != label)
                    {
                        if(Labels[node] == 1)
                        {
                            labeling[node] = label;
                            nodesChanged++;
                        }
                    }
                }

                if(verbose)
                {
                    energy->report();
                    std::cout << "  LAB " << label << ":\t" << lastEnergy << " -> " << newEnergy << " / " << ratioUnlabeled * 100 << "% UNL / "
                    << nodesChanged / static_cast<double>(energy->getNumNodes()) * 100 << "% CHN / IMP: " << improveCounter << std::endl;
                    lastEnergy = newEnergy;
                }
            }
        }
    }
#ifndef PRINT_ENERGY
    lastEnergy = energy->evaluateTotalCostSum();
#endif
    return lastEnergy;
}

} //namespace newmeshreg
