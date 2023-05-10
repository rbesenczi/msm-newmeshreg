#ifndef NEWMESHREG_FUSION_H
#define NEWMESHREG_FUSION_H

#ifdef HAS_HOCR
#include "ELC/ELC.h"
#endif

#include "DiscreteCostFunction.h"
#include "FastPD.h"

#define NUM_SWEEPS 2

typedef	double REAL;
struct UnaryData	{ REAL buffer[2]; };
struct TripletData	{ REAL buffer[8]; };

namespace ELCReduce { template<typename T> class PBF; }

namespace newmeshreg {

class DummyCostFunction: public DiscreteCostFunction {

public:
    DummyCostFunction(){ m_num_labels = 2; }

    void setUnaryCost(int node, double cost0, double cost1){
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

    void set_parameters(myparam&){}

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

enum Reduction { ELC_HOCR, ELC_APPROX, HOCR };

class Fusion {
public:

    template<typename OPTIMIZER>
    static void reduce_and_convert(ELCReduce::PBF<REAL>& pbf, OPTIMIZER& MODEL, Reduction mode) {

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

    static double optimize(const std::shared_ptr<DiscreteModel>& energy, Reduction reductionMode, bool verbose, int numthreads) {

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
                sumlabeldiff = 0.0;
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
};

} //namespace newmeshreg

#endif //NEWMESHREG_FUSION_H
