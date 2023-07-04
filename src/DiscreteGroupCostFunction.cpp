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
    it=p.find("quartet");_quadcost=boost::get<bool>(it->second);
}

void DiscreteGroupCostFunction::get_spacings() {

}

void DiscreteGroupCostFunction::initialize(int numNodes, int numLabels, int numPairs, int numTriplets,
                                           int numQuartets) {

}

void DiscreteGroupCostFunction::define_template_patches() {

}

void DiscreteGroupCostFunction::resample_to_template() {

}

double DiscreteGroupCostFunction::computeQuartetCost(int quartet, int labelA, int labelB, int labelC, int labelD) {
    return 0;
}

double DiscreteGroupCostFunction::computeTripletCost(int triplet, int labelA, int labelB, int labelC) {
    return 0;
}

double DiscreteGroupCostFunction::computePairwiseCost(int pair, int labelA, int labelB) {
    return 0;
}

void DiscreteGroupCostFunction::resample_patches() {

}

void DiscreteGroupCostFunction::resampler_worker_function(int begin, int end, const std::vector<bool>& INrange) {

}

std::map<int,float> DiscreteGroupCostFunction::resample_onto_template(int node, int ID, const newresampler::Point& newCP,
                                                                      const std::vector<int>& PTS) {

}

} //namespace newmeshreg
