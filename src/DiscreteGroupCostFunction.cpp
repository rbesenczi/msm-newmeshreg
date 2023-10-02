#include "DiscreteGroupCostFunction.h"

namespace newmeshreg {

void DiscreteGroupCostFunction::initialize(int numNodes, int numLabels, int numPairs, int numTriplets) {
    DiscreteCostFunction::initialize(numNodes, numLabels, numPairs, numTriplets);
    _sourceinrange.clear();
    _sourceinrange.resize(VERTICES_PER_SUBJ * num_subjects * numLabels);
    get_source_data();
}

double DiscreteGroupCostFunction::computeTripletCost(int triplet, int labelA, int labelB, int labelC) {

    int meshID = std::floor(triplet/(double)TRIPLETS_PER_SUBJ);
    int m0 = _triplets[3*triplet]-meshID*VERTICES_PER_SUBJ;
    int m1 = _triplets[3*triplet+1]-meshID*VERTICES_PER_SUBJ;
    int m2 = _triplets[3*triplet+2]-meshID*VERTICES_PER_SUBJ;

    newresampler::Triangle TRI((*ROTATIONS)[_triplets[3*triplet  ]] * _labels[labelA],
                               (*ROTATIONS)[_triplets[3*triplet+1]] * _labels[labelB],
                               (*ROTATIONS)[_triplets[3*triplet+2]] * _labels[labelC],0);
    newresampler::Triangle TRI_noDEF(_CONTROLMESHES[meshID].get_coord(m0),
                                     _CONTROLMESHES[meshID].get_coord(m1),
                                     _CONTROLMESHES[meshID].get_coord(m2),0);

    if((TRI.normal() | TRI_noDEF.normal()) < 0) { return 1e7; }

    newresampler::Triangle TRI_ORIG(_ORIG.get_coord(m0),
                                    _ORIG.get_coord(m1),
                                    _ORIG.get_coord(m2),0);

    return _reglambda * MISCMATHS::pow(calculate_triangular_strain(TRI_ORIG,TRI,_mu,_kappa),_rexp);
}

void DiscreteGroupCostFunction::get_source_data() {
    for (int subject = 0; subject < num_subjects; subject++)
        #pragma omp parallel for num_threads(_threads)
        for (int vertex = 0; vertex < VERTICES_PER_SUBJ; vertex++)
            for (int label = 0; label < m_num_labels; label++) {
                newresampler::Point LP = (*ROTATIONS)[subject * VERTICES_PER_SUBJ + vertex] * _labels[label];
                for (int datapoint = 0; datapoint < _DATAMESHES[subject].nvertices(); datapoint++)
                    if (((2 * RAD * asin((LP - _DATAMESHES[subject].get_coord(datapoint)).norm() / (2 * RAD))) <
                         _controlptrange * SPACINGS[subject](vertex + 1)))
                _sourceinrange.at(subject * VERTICES_PER_SUBJ * m_num_labels + vertex * m_num_labels + label).push_back(datapoint);
            }
}

double DiscreteGroupCostFunction::computePairwiseCost(int pair, int labelA, int labelB) {

    int nodeA = _pairs[2*pair];
    int nodeB = _pairs[2*pair+1];

    int meshAID = std::floor((double)nodeA/VERTICES_PER_SUBJ);
    int meshBID = std::floor((double)nodeB/VERTICES_PER_SUBJ);

    int nodeAID = nodeA - meshAID*VERTICES_PER_SUBJ;
    int nodeBID = nodeB - meshBID*VERTICES_PER_SUBJ;

    std::vector<int> patchA = _sourceinrange[meshAID*VERTICES_PER_SUBJ*m_num_labels + nodeAID*m_num_labels + labelA];
    std::vector<int> patchB = _sourceinrange[meshBID*VERTICES_PER_SUBJ*m_num_labels + nodeBID*m_num_labels + labelB];

    std::vector<double> patch_data_A(patchA.size(),0.0),
                        patch_data_B(patchB.size(),0.0);

    for(int i = 0; i < patchA.size(); ++i)
        patch_data_A[i] = FEAT->get_data_val(1,patchA[i] + 1,meshAID);
    for(int i = 0; i < patchB.size(); ++i)
        patch_data_B[i] = FEAT->get_data_val(1,patchB[i] + 1,meshBID);

    if(patch_data_A.size() < patch_data_B.size())
        patch_data_B.resize(patch_data_A.size());
    else if (patch_data_A.size() > patch_data_B.size())
        patch_data_A.resize(patch_data_B.size());

    return sim.get_sim_for_min(patch_data_A, patch_data_B);
}

} //namespace newmeshreg
