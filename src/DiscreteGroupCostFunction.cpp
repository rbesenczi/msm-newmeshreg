#include "DiscreteGroupCostFunction.h"

namespace newmeshreg {

void DiscreteGroupCostFunction::initialize(int numNodes, int numLabels, int numPairs, int numTriplets) {
    DiscreteCostFunction::initialize(numNodes, numLabels, numPairs, numTriplets);
    source_in_range_data.clear();
    source_in_range_data.resize(VERTICES_PER_SUBJ * num_subjects * numLabels);
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
                for (int label = 0; label < m_num_labels; label++)
                {
                    std::map<int, double> patchdata;
                    const newresampler::Point CP = _CONTROLMESHES[subject].get_coord(vertex);
                    for (int datapoint = 0; datapoint < _DATAMESHES[subject].nvertices(); datapoint++)
                        if (((2 * RAD * asin((CP - _DATAMESHES[subject].get_coord(datapoint)).norm() / (2 * RAD))) <
                             _controlptrange * SPACINGS[subject](vertex + 1)))
                        {
                            newresampler::Point tmp =
                                    estimate_rotation_matrix(CP, (*ROTATIONS)[subject * VERTICES_PER_SUBJ + vertex] * _labels[label])
                                    * _DATAMESHES[subject].get_coord(datapoint);

                            newresampler::Triangle closest_triangle = targettree->get_closest_triangle(tmp);
                            int closest_node_id = targettree->get_closest_vertex_ID(tmp);

                            newresampler::Point
                                    v0 = closest_triangle.get_vertex_coord(0),
                                    v1 = closest_triangle.get_vertex_coord(1),
                                    v2 = closest_triangle.get_vertex_coord(2);
                                int n0 = closest_triangle.get_vertex_no(0),
                                    n1 = closest_triangle.get_vertex_no(1),
                                    n2 = closest_triangle.get_vertex_no(2);

                            newresampler::project_point(v0,v1,v2,tmp);

                            /* TODO need something like this...
                             * SP=_DATAMESHES[ID].get_coord(sourceind-1);
	                            projectPoint(CP0,CP1,CP2,SP);
	                            tmp= barycentric(CP0,CP1,CP2,SP,newCP0,newCP1,newCP2);
                             */

                            patchdata[closest_node_id] = newresampler::barycentric_weight(v0, v1, v2, tmp,
                                                                                FEAT->get_data_val(1, datapoint + 1,subject),
                                                                                FEAT->get_data_val(1, n1 + 1),
                                                                                FEAT->get_data_val(1, n2 + 1));
                            source_in_range_data[subject * VERTICES_PER_SUBJ * m_num_labels + vertex * m_num_labels +
                                                 label].push_back(FEAT->get_data_val(1, datapoint + 1, subject));
                        }
                }

}

double DiscreteGroupCostFunction::computePairwiseCost(int pair, int labelA, int labelB) {

    int subject_A = std::floor((double)_pairs[2*pair  ] / VERTICES_PER_SUBJ);
    int subject_B = std::floor((double)_pairs[2*pair+1] / VERTICES_PER_SUBJ);

    int node_A = _pairs[2*pair  ] - subject_A * VERTICES_PER_SUBJ;
    int node_B = _pairs[2*pair+1] - subject_B * VERTICES_PER_SUBJ;

    std::vector<double> patchA = source_in_range_data[subject_A * VERTICES_PER_SUBJ * m_num_labels + node_A * m_num_labels + labelA];
    std::vector<double> patchB = source_in_range_data[subject_B * VERTICES_PER_SUBJ * m_num_labels + node_B * m_num_labels + labelB];

    // patch sizes must be equal for Pearson's correlation
    if(patchA.size() != patchB.size()) {
        if (patchA.size() < patchB.size()) patchB.resize(patchA.size());
        if (patchA.size() > patchB.size()) patchA.resize(patchB.size());
    }

    return sim.get_sim_for_min(patchA, patchB);
}

} //namespace newmeshreg
