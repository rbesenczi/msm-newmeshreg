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

double Fusion::optimize(const std::shared_ptr<DiscreteModel>& energy, Reduction reductionMode, bool verbose) {

    const int numNodes    = energy->getNumNodes();
    const int numPairs    = energy->getNumPairs();
    const int numTriplets = energy->getNumTriplets();
    const int *pairs    = energy->getPairs();
    const int *triplets = energy->getTriplets();

#ifdef HAS_QPBO
    if(verbose) cout << " has QPBO " << endl;
    QPBO<REAL> qpbo(numGraphNodes,numGraphEdges) ;
#endif
    std::shared_ptr<DiscreteModelDummy> FPDMODEL = std::make_shared<DiscreteModelDummy>();

    const int numSweeps = NUM_SWEEPS;
    const int numLabels = energy->getNumLabels();
    int *labeling = energy->getLabeling();

    double initEnergy = energy->evaluateTotalCostSum();

    double lastEnergy = initEnergy;
    double unlabeledAvg = 0.0;
    double sumlabeldiff = 0.0;

    for(int sweep = 0; sweep < numSweeps; ++sweep)
    {
        for(int label = 0; label < numLabels; ++label)
        {
            sumlabeldiff = 0;
            ELCReduce::PBF<REAL> pbf;
            int improveCounter = 0;
            double ratioUnlabeled = 0;
            int nodesChanged = 0;
            std::vector<UnaryData> unary_data(numNodes);

            for(int node = 0; node < numNodes; ++node)
            {
                unary_data[node].buffer[0] = energy->computeUnaryCost(node,labeling[node]);	//0
                unary_data[node].buffer[1] = energy->computeUnaryCost(node,label);				//1
                sumlabeldiff += abs(label-labeling[node]);
            }

            if(sumlabeldiff > 0)
            {
                for (int node = 0; node < numNodes; ++node)
                    pbf.AddUnaryTerm(node, unary_data[node].buffer[0], unary_data[node].buffer[1]);

                std::vector<PairData> pair_data(numPairs);

                #pragma omp parallel for
                for(int pair = 0; pair < numPairs; ++pair)
                {
                    const int nodeA = pairs[pair*2];
                    const int nodeB = pairs[pair*2+1];

                    pair_data[pair].buffer[0] = energy->computePairwiseCost(pair,labeling[nodeA],labeling[nodeB]);	//00
                    pair_data[pair].buffer[1] = energy->computePairwiseCost(pair,labeling[nodeA],label);			//01
                    pair_data[pair].buffer[2] = energy->computePairwiseCost(pair,label,labeling[nodeB]);			//10
                    pair_data[pair].buffer[3] = energy->computePairwiseCost(pair,label,label);
                }

                for(int pair = 0; pair < numPairs; ++pair)
                {
                    const int nodeA = pairs[pair*2];
                    const int nodeB = pairs[pair*2+1];
                    int node_ids[2] = { nodeA, nodeB};

                    pbf.AddPairwiseTerm(node_ids[0], node_ids[1], pair_data[pair].buffer[0], pair_data[pair].buffer[1], pair_data[pair].buffer[2], pair_data[pair].buffer[3]);
                }

                std::vector<TripletData> triplet_data(numTriplets);

                #pragma omp parallel for
                for (int triplet = 0; triplet < numTriplets; ++triplet)
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

                for (int triplet = 0; triplet < numTriplets; ++triplet)
                {
                    const int nodeA = triplets[triplet*3];
                    const int nodeB = triplets[triplet*3+1];
                    const int nodeC = triplets[triplet*3+2];
                    int node_ids[3] = { nodeA, nodeB, nodeC };
                    pbf.AddHigherTerm(3, node_ids, triplet_data[triplet].buffer);
                }

                double newEnergy;
#ifdef HAS_QPBO
                qpbo.Reset();
                reduce_and_convert(pbf, qpbo, reductionMode);
                qpbo.MergeParallelEdges();
                qpbo.Solve();

                for(int node = 0; node < numNodes; ++node)
                {
                    if(labeling[node] != label)
                    {
                        if(qpbo.GetLabel(node) < 0) unlabeledNodes++;
                    }
                }

                //TRY QPBO-I to improve the solution
                if(unlabeledNodes > 0)
                {
                    srand ( static_cast<unsigned int>(time(NULL)) );

                    int numTrials = MAX_IMPROVEMENTS;
                    const double ratioThresh = 0.3;
                    ratioUnlabeled = static_cast<double>(unlabeledNodes) / static_cast<double>(numNodes);
                    if (ratioUnlabeled < ratioThresh) numTrials = static_cast<int>(0.5+ratioUnlabeled*numTrials/ratioThresh);

                    if(MAX_IMPROVEMENTS > 0)
                    {
                        for(int i = 0; i < numTrials; ++i)
                        {
                            if(qpbo.Improve()) improveCounter++;
                        }
                    }

                    if(ratioUnlabeled > unlabeledMax) unlabeledMax = ratioUnlabeled;
                    unlabeledAvg += ratioUnlabeled;
                }

                for(int node = 0; node < numNodes; ++node)
                {
                    if(labeling[node] != label)
                    {
                        if(qpbo.GetLabel(node) == 1)
                        {
                            labeling[node] = label;
                            nodesChanged++;
                        }
                    }
                }

                if(verbose)
                {
                    //energy->applytestLabeling(labeling,label);
                    newEnergy = energy->evaluateTotalCostSum();
                }
#else
                FPDMODEL->reset();
                reduce_and_convert(pbf,*FPDMODEL, reductionMode);
                FPDMODEL->initialise();
                int *Labels = FPDMODEL->getLabeling();
                FPD::FastPD opt(FPDMODEL, 100 );
                newEnergy = opt.run();
                opt.getLabeling(Labels);

                for(int node = 0; node < numNodes; ++node)
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
#endif
                if(verbose)
                {
                    energy->report();
                    std::cout << "  LAB " << label << ":\t" << lastEnergy << " -> " << newEnergy << " / " << ratioUnlabeled * 100 << "% UNL / "
                    << nodesChanged / static_cast<double>(numNodes) * 100 << "% CHN / IMP: " << improveCounter << std::endl;
                    lastEnergy = newEnergy;
                }
            }
        }
    }


    unlabeledAvg /= static_cast<double>(numLabels*numSweeps);
#ifndef PRINT_ENERGY
    lastEnergy = energy->evaluateTotalCostSum();
#endif
    return lastEnergy;
}

} //namespace newmeshreg
