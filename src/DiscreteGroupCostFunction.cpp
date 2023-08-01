#include "DiscreteGroupCostFunction.h"

namespace newmeshreg {

void DiscreteGroupCostFunction::set_parameters(myparam& p) {
    myparam::iterator it;
    it=p.find("lambda");_reglambda=boost::get<float>(it->second);
    it=p.find("lambda_pairs");_lambdapairs=boost::get<float>(it->second);
    it=p.find("set_lambda_pairs");_setpairs=boost::get<bool>(it->second);
    it=p.find("simmeasure");_simmeasure=boost::get<int>(it->second); sim.set_simval(_simmeasure);
    it=p.find("verbosity");_verbosity=boost::get<bool>(it->second);
    it=p.find("shearmodulus");_mu=boost::get<float>(it->second);
    it=p.find("bulkmodulus");_kappa=boost::get<float>(it->second);
    it=p.find("sigma_in");_sigma=boost::get<float>(it->second);
    it=p.find("numthreads"); _threads=boost::get<int>(it->second);
}

void DiscreteGroupCostFunction::initialize(int numNodes, int numLabels, int numPairs, int numTriplets) {
    DiscreteCostFunction::initialize(numNodes, numLabels, numPairs, numTriplets);
    resample_to_template();
    get_spacings();
    _sourceinrange.clear();
    _sourceinrange.resize(VERTICES_PER_SUBJ * num_subjects);
}

void DiscreteGroupCostFunction::get_spacings() {

    NEWMAT::ColumnVector vMAXmvd(_CONTROLMESHES[0].nvertices());
    SPACINGS.clear();

    for(int n = 0; n < num_subjects; n++)
    {
        vMAXmvd = 0;
        for (int k = 0; k < _CONTROLMESHES[n].nvertices(); k++)
        {
            newresampler::Point CP = _CONTROLMESHES[n].get_coord(k);
            for (auto it = _CONTROLMESHES[n].nbegin(k); it != _CONTROLMESHES[n].nend(k); it++)
            {
                double dist = 2 * RAD * asin((CP -_CONTROLMESHES[n].get_coord(*it)).norm() / (2*RAD));
                if(dist > vMAXmvd(k+1)) vMAXmvd(k+1) = dist;
            }
        }
        SPACINGS.push_back(vMAXmvd);
    }
}

void DiscreteGroupCostFunction::resample_to_template() {

    RESAMPLEDDATA.clear();
    RESAMPLEDDATA.resize(num_subjects);

    for(int n = 0; n < num_subjects; n++)
    {
        _DATAMESHES[n].set_pvalues(FEAT->get_data_matrix(n));
        RESAMPLEDDATA[n] = newresampler::metric_resample(_DATAMESHES[n], _TEMPLATE).get_pvalues();
    }
}

double DiscreteGroupCostFunction::computeTripletCost(int triplet, int labelA, int labelB, int labelC) {

    int meshID = floor(triplet/TRIPLETS_PER_SUBJ);
    int m0 = _triplets[3*triplet]-meshID*VERTICES_PER_SUBJ;
    int m1 = _triplets[3*triplet+1]-meshID*VERTICES_PER_SUBJ;
    int m2 = _triplets[3*triplet+2]-meshID*VERTICES_PER_SUBJ;

    newresampler::Triangle TRI((*ROTATIONS)[_triplets[3*triplet  ]] * _labels[labelA],
                               (*ROTATIONS)[_triplets[3*triplet+1]] * _labels[labelB],
                               (*ROTATIONS)[_triplets[3*triplet+2]] * _labels[labelC],0);
    newresampler::Triangle TRI_ORIG(_ORIG.get_coord(m0),
                                    _ORIG.get_coord(m1),
                                    _ORIG.get_coord(m2),0);
    newresampler::Triangle TRI_noDEF(_CONTROLMESHES[meshID].get_coord(m0),
                                     _CONTROLMESHES[meshID].get_coord(m1),
                                     _CONTROLMESHES[meshID].get_coord(m2),0);

    if((TRI.normal() | TRI_noDEF.normal()) < 0) { return 1e7; }

    return _reglambda * MISCMATHS::pow(calculate_triangular_strain(TRI_ORIG,TRI,_mu,_kappa),_rexp);
}

void DiscreteGroupCostFunction::get_source_data() {
    for (int subject = 0; subject < num_subjects; ++subject)
        for (int k = 0; k < _CONTROLMESHES[subject].nvertices(); k++)
            for (int i = 0; i < _DATAMESHES[subject].nvertices(); i++)
                if (((2 * RAD * asin((_CONTROLMESHES[subject].get_coord(k)-_DATAMESHES[subject].get_coord(i)).norm() / (2*RAD)))
                            < _controlptrange * SPACINGS[subject](k+1)))
                    _sourceinrange[subject*VERTICES_PER_SUBJ+k].push_back(i);
}

double DiscreteGroupCostFunction::computePairwiseCost(int pair, int labelA, int labelB) {

    int nodeA = _pairs[2*pair];
    int nodeB = _pairs[2*pair+1];

    int meshAID = std::floor((double)nodeA/VERTICES_PER_SUBJ);
    int meshBID = std::floor((double)nodeB/VERTICES_PER_SUBJ);

    int nodeAID = nodeA - meshAID*VERTICES_PER_SUBJ;
    int nodeBID = nodeB - meshBID*VERTICES_PER_SUBJ;

    auto rotA = newresampler::estimate_rotation_matrix(_CONTROLMESHES[meshAID].get_coord(nodeAID), (*ROTATIONS)[nodeA]*_labels[labelA]);
    auto rotB = newresampler::estimate_rotation_matrix(_CONTROLMESHES[meshBID].get_coord(nodeBID), (*ROTATIONS)[nodeB]*_labels[labelB]);

    return sim.get_sim_for_min(get_patch_data(nodeA, rotA),
                                get_patch_data(nodeB, rotB));
}

std::vector<double> DiscreteGroupCostFunction::get_patch_data(int node, const NEWMAT::Matrix& rot) {

    std::vector<double> data(_sourceinrange[node].size(), 0.0);
    int subject = std::floor((double)node/VERTICES_PER_SUBJ);

    #pragma omp parallel for num_threads(_threads)
    for(unsigned int i = 0; i < _sourceinrange[node].size(); i++)
    {
        newresampler::Point tmp = rot * _DATAMESHES[subject].get_coord(_sourceinrange[node][i]);

        newresampler::Triangle closest_triangle = datameshtrees[subject]->get_closest_triangle(tmp);

        newresampler::Point v0 = closest_triangle.get_vertex_coord(0),
                            v1 = closest_triangle.get_vertex_coord(1),
                            v2 = closest_triangle.get_vertex_coord(2);
                        int n0 = closest_triangle.get_vertex_no(0),
                            n1 = closest_triangle.get_vertex_no(1),
                            n2 = closest_triangle.get_vertex_no(2);

        data[i] = newresampler::barycentric_weight(v0, v1, v2, tmp,
                                                    FEAT->get_data_val(1, n0+1, subject),
                                                    FEAT->get_data_val(1, n1+1, subject),
                                                    FEAT->get_data_val(1, n2+1, subject));
    }

    return data;
}

} //namespace newmeshreg
