#ifndef NEWMESHREG_FUSION_H
#define NEWMESHREG_FUSION_H

#include "ELC/ELC.h"
#include "../src/DiscreteCostFunction.h"
#include "FastPD/FastPD.h"

#define NUM_SWEEPS 2

namespace newmeshreg {

struct UnaryData    { double buffer[2]; };
struct PairData     { double buffer[4]; };
struct TripletData  { double buffer[8]; };
struct QuartetData  { double buffer[16]; };

class DummyCostFunction: public DiscreteCostFunction {

public:
    DummyCostFunction() { m_num_labels = 2; }

    void setUnaryCost(int node, double cost0, double cost1) {
        unaryenergies.insert(std::pair<int, std::vector<double>>(node, std::vector<double>()));
        unaryenergies[node].push_back(cost0);
        unaryenergies[node].push_back(cost1);
    }

    void setPairwiseCost(int ind, double E00, double E01, double E10, double E11) {
        pairenergies.insert(std::pair<int, std::vector<double>>(ind, std::vector<double>()));
        pairenergies[ind].push_back(E00);
        pairenergies[ind].push_back(E01);
        pairenergies[ind].push_back(E10);
        pairenergies[ind].push_back(E11);
    }

    double computePairwiseCost(int pair, int labelA, int labelB) override {
        if(labelA == 0 && labelB == 0)
            return pairenergies[pair][0];
        else if (labelA==0 && labelB==1)
            return pairenergies[pair][1];
        else if (labelA==1 && labelB==0)
            return pairenergies[pair][2];
        else
            return pairenergies[pair][3];
    }

    void convertenergies(int numNodes,int numPairs, int numLabels) {

        m_num_nodes = numNodes;
        m_num_labels = numLabels;
        m_num_pairs = numPairs;

        delete[] unarycosts;
        unarycosts = new double[numNodes*numLabels];

        for (int i = 0; i < numLabels; i++)
            for (int j = 0; j < numNodes; j++)
                unarycosts[i * numNodes + j] = unaryenergies[j][i];
    }

    void reset() override {
        unaryenergies.clear();
        pairenergies.clear();
    }

protected:
    std::map<int,std::vector<double>> unaryenergies; // maps of nodes  xlabels x vals
    std::map<int,std::vector<double>> pairenergies;
};

class DiscreteModelDummy : public DiscreteModel {

public:
    DiscreteModelDummy() {
        costfct = std::make_shared<DummyCostFunction>();
        m_num_pairs = 0; m_num_nodes = 0; m_num_labels = 2;
    }

    // create dummy costfunction to be used with dummy model (will just save output of conversion)
    std::shared_ptr<DiscreteCostFunction> getCostFunction() override { return costfct; }

    //---ELC conversion functions---//
    void AddNode(int num){ m_num_nodes = num; }

    // Adds unary term Ei(x_i) to the energy function with cost values Ei(0)=E0, Ei(1)=E1.
    void AddUnaryTerm(int node, double E0, double E1) { costfct->setUnaryCost(node,E0, E1); }

    // adds pairwise term for binary costs with label combinations 00,01,10,11
    void AddPairwiseTerm(int node1, int node2, double E00, double E01, double E10, double E11) {
        pairIDs.insert(std::pair<int,std::vector<int>>(m_num_pairs, std::vector<int>()));
        pairIDs[m_num_pairs].push_back(node1);
        pairIDs[m_num_pairs].push_back(node2);
        costfct->setPairwiseCost(m_num_pairs,E00,E01,E10,E11);
        m_num_pairs++;
    }

    // FastPD conversion functions
    void initialise(){
        initLabeling();
        costfct->convertenergies(m_num_nodes,m_num_pairs,2);
        pairs = new int[m_num_pairs * 2];
        for(int i = 0; i < m_num_pairs; i++)
        {
            pairs[2 * i] = pairIDs[i][0];
            pairs[2 * i + 1] = pairIDs[i][1];
        }
    }

    void reset() {
        pairIDs.clear();
        m_num_pairs = 0;
        m_num_nodes = 0;
        m_num_labels = 2;
        costfct->reset();
    }

protected:
    std::map<int,std::vector<int>> pairIDs;
    std::shared_ptr<DummyCostFunction> costfct;
};

class Fusion {
public:
    static double optimize(const std::shared_ptr<DiscreteModel>& energy, bool verbose, int numthreads) {

        const int *pairs    = energy->getPairs();
        const int* triplets = energy->getTriplets();
        const int *quartets = energy->getQuartets();

        int *labeling = energy->getLabeling();
        std::shared_ptr<DiscreteModelDummy> FPDMODEL = std::make_shared<DiscreteModelDummy>();

        double initEnergy = energy->evaluateTotalCostSum();
        double lastEnergy = initEnergy;
        double sumlabeldiff = 0.0;

        for(int sweep = 0; sweep < NUM_SWEEPS; ++sweep)
        {
            for(int label = 0; label < energy->getNumLabels(); ++label)
            {
                sumlabeldiff = 0.0;
                ELCReduce::PBF<double> pbf;
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
                    for(int node = 0; node < energy->getNumNodes(); ++node)
                        pbf.AddUnaryTerm(node, unary_data[node].buffer[0], unary_data[node].buffer[1]);

                    std::vector<PairData> pair_data(energy->getNumPairs());

                    #pragma omp parallel for num_threads(numthreads)
                    for(int pair = 0; pair < energy->getNumPairs(); ++pair)
                    {
                        const int nodeA = pairs[pair*2];
                        const int nodeB = pairs[pair*2+1];

                        pair_data[pair].buffer[0] = energy->computePairwiseCost(pair,labeling[nodeA],labeling[nodeB]);    //00
                        pair_data[pair].buffer[1] = energy->computePairwiseCost(pair,labeling[nodeA],label);          //01
                        pair_data[pair].buffer[2] = energy->computePairwiseCost(pair,label,labeling[nodeB]);          //10
                        pair_data[pair].buffer[3] = energy->computePairwiseCost(pair,label,label); //11
                    }

                    for(int pair = 0; pair < energy->getNumPairs(); ++pair)
                    {
                        const int nodeA = pairs[pair*2];
                        const int nodeB = pairs[pair*2+1];
                        int node_ids[2] = { nodeA, nodeB };
                        pbf.AddPairwiseTerm(node_ids[0], node_ids[1], pair_data[pair].buffer[0], pair_data[pair].buffer[1], pair_data[pair].buffer[2], pair_data[pair].buffer[3]);
                    }

                    std::vector<TripletData> triplet_data(energy->getNumTriplets());

                    #pragma omp parallel for num_threads(numthreads)
                    for(int triplet = 0; triplet < energy->getNumTriplets(); ++triplet)
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

                    std::vector<QuartetData> quartet_data(energy->getNumQuartets());

                    #pragma omp parallel for num_threads(numthreads)
                    for (int quartet = 0; quartet < energy->getNumQuartets(); ++quartet)
                    {
                        const int nodeA = quartets[quartet*4];
                        const int nodeB = quartets[quartet*4+1];
                        const int nodeC = quartets[quartet*4+2];
                        const int nodeD = quartets[quartet*4+3];
                        if(!(label==labeling[nodeA] && label==labeling[nodeB] && label==labeling[nodeC] && label==labeling[nodeD]))
                        {
                            quartet_data[quartet].buffer[0]  = energy->computeQuartetCost(quartet,labeling[nodeA],labeling[nodeB],labeling[nodeC],labeling[nodeD]);
                            quartet_data[quartet].buffer[1]  = energy->computeQuartetCost(quartet,labeling[nodeA],labeling[nodeB],labeling[nodeC],label);
                            quartet_data[quartet].buffer[2]  = energy->computeQuartetCost(quartet,labeling[nodeA],labeling[nodeB],label,labeling[nodeD]);
                            quartet_data[quartet].buffer[3]  = energy->computeQuartetCost(quartet,labeling[nodeA],labeling[nodeB],label,label);
                            quartet_data[quartet].buffer[4]  = energy->computeQuartetCost(quartet,labeling[nodeA],label,labeling[nodeC],labeling[nodeD]);
                            quartet_data[quartet].buffer[5]  = energy->computeQuartetCost(quartet,labeling[nodeA],label,labeling[nodeC],label);
                            quartet_data[quartet].buffer[6]  = energy->computeQuartetCost(quartet,labeling[nodeA],label,label,labeling[nodeD]);
                            quartet_data[quartet].buffer[7]  = energy->computeQuartetCost(quartet,labeling[nodeA],label,label,label);
                            quartet_data[quartet].buffer[8]  = energy->computeQuartetCost(quartet,label,labeling[nodeB],labeling[nodeC],labeling[nodeD]);
                            quartet_data[quartet].buffer[9]  = energy->computeQuartetCost(quartet,label,labeling[nodeB],labeling[nodeC],label);
                            quartet_data[quartet].buffer[10] = energy->computeQuartetCost(quartet,label,labeling[nodeB],label,labeling[nodeD]);
                            quartet_data[quartet].buffer[11] = energy->computeQuartetCost(quartet,label,labeling[nodeB],label,label);
                            quartet_data[quartet].buffer[12] = energy->computeQuartetCost(quartet,label,label,labeling[nodeC],labeling[nodeD]);
                            quartet_data[quartet].buffer[13] = energy->computeQuartetCost(quartet,label,label,labeling[nodeC],label);
                            quartet_data[quartet].buffer[14] = energy->computeQuartetCost(quartet,label,label,label,labeling[nodeD]);
                            quartet_data[quartet].buffer[15] = energy->computeQuartetCost(quartet,label,label,label,label);
                        }
                    }

                    for (int quartet = 0; quartet < energy->getNumQuartets(); ++quartet)
                    {
                        const int nodeA = quartets[quartet*4];
                        const int nodeB = quartets[quartet*4+1];
                        const int nodeC = quartets[quartet*4+2];
                        const int nodeD = quartets[quartet*4+3];
                        if(!(label==labeling[nodeA] && label==labeling[nodeB] && label==labeling[nodeC] && label==labeling[nodeD]))
                        {
                            int node_ids[4] = { nodeA, nodeB, nodeC, nodeD };
                            pbf.AddHigherTerm(4, node_ids, quartet_data[quartet].buffer);
                        }
                    }

                    FPDMODEL->reset();

                    ELCReduce::PBF<double> qpbf;
                    pbf.toQuadratic(qpbf, pbf.maxID()+1); // Reduce to Quadratic pseudo-Boolean function using HOCR.
                    qpbf.convert(*FPDMODEL, qpbf.maxID()+1);
                    pbf.clear();
                    qpbf.clear();

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
        lastEnergy = energy->evaluateTotalCostSum();

        return lastEnergy;
    }
};

} //namespace newmeshreg

#endif //NEWMESHREG_FUSION_H
