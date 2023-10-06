#include "DiscreteGroupCostFunction.h"

namespace newmeshreg {

void DiscreteGroupCostFunction::initialize(int numNodes, int numLabels, int numPairs, int numTriplets) {
    DiscreteCostFunction::initialize(numNodes, numLabels, numPairs, numTriplets);
    patch_data.clear();
    patch_data.resize(VERTICES_PER_SUBJ * num_subjects * numLabels);
    get_patch_data();
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

void DiscreteGroupCostFunction::get_patch_data() {
    for (int subject = 0; subject < num_subjects; subject++)
        #pragma omp parallel for num_threads(_threads)
        for (int vertex = 0; vertex < VERTICES_PER_SUBJ; vertex++)
        {
            const newresampler::Point CP = _CONTROLMESHES[subject].get_coord(vertex);
            for (int label = 0; label < m_num_labels; label++)
            {
                std::map<int, double> patchdata;
                for (int datapoint = 0; datapoint < _DATAMESHES[subject].nvertices(); datapoint++)
                    if (((2 * RAD * asin((CP - _DATAMESHES[subject].get_coord(datapoint)).norm() / (2 * RAD))) <
                         _controlptrange * SPACINGS[subject](vertex + 1))) {
                        newresampler::Point tmp =
                                estimate_rotation_matrix(CP, (*ROTATIONS)[subject * VERTICES_PER_SUBJ + vertex] * _labels[label])
                                                * _DATAMESHES[subject].get_coord(datapoint);
                        patchdata[targettree->get_closest_vertex_ID(tmp)] = FEAT->get_data_val(1, datapoint + 1, subject);
                    }
                patch_data[subject * VERTICES_PER_SUBJ * m_num_labels + vertex * m_num_labels + label] = patchdata;
            }
        }
}

double DiscreteGroupCostFunction::computePairwiseCost(int pair, int labelA, int labelB) {

    std::vector<double> patch_data_A, patch_data_B;

    int subject_A = std::floor((double)_pairs[2*pair  ] / VERTICES_PER_SUBJ);
    int subject_B = std::floor((double)_pairs[2*pair+1] / VERTICES_PER_SUBJ);

    int node_A = _pairs[2*pair  ] - subject_A * VERTICES_PER_SUBJ;
    int node_B = _pairs[2*pair+1] - subject_B * VERTICES_PER_SUBJ;

    std::map<int,double> patchA = patch_data[subject_A * VERTICES_PER_SUBJ * m_num_labels + node_A * m_num_labels + labelA];
    std::map<int,double> patchB = patch_data[subject_B * VERTICES_PER_SUBJ * m_num_labels + node_B * m_num_labels + labelB];

    for (const auto &e: patchA)
        if (patchB.find(e.first) != patchB.end())
            patch_data_A.push_back(e.second);

    for (const auto &e: patchB)
        if (patchA.find(e.first) != patchA.end())
            patch_data_B.push_back(e.second);

    // patch sizes must be equal for Pearson's correlation
    if(patch_data_A.size() != patch_data_B.size()) {
        std::cout << "Is it?" << std::endl;
        if (patch_data_A.size() < patch_data_B.size()) patch_data_B.resize(patch_data_A.size());
        if (patch_data_A.size() > patch_data_B.size()) patch_data_A.resize(patch_data_B.size());
    }

    return sim.get_sim_for_min(patch_data_A, patch_data_B);
}

} //namespace newmeshreg
