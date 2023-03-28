#include "DiscreteCostFunction.h"

using namespace std;

namespace newmeshreg {

//================================BASE CLASS===========================================================================//
void DiscreteCostFunction::initialize(int numNodes, int numLabels, int numPairs, int numTriplets) {

    if (m_num_nodes != numNodes || m_num_labels != numLabels)
    {
        delete[] unarycosts;
        unarycosts = new double[numNodes * numLabels];
    }

    if (m_num_pairs != numPairs || m_num_labels != numLabels)
    {
        delete[] paircosts;
        paircosts = new double[numPairs * numLabels * numLabels];
    }

    m_num_nodes = numNodes;
    m_num_labels = numLabels;
    m_num_pairs = numPairs;
    m_num_triplets = numTriplets;

    std::fill(unarycosts,unarycosts+m_num_labels*m_num_nodes,0.0f);
    std::fill(paircosts,paircosts+m_num_labels*m_num_labels*m_num_pairs,0.0f);
}

void DiscreteCostFunction::reset() {
    if(unarycosts) std::fill(unarycosts,unarycosts+m_num_labels*m_num_nodes,0.0f);
    if(paircosts) std::fill(paircosts,paircosts+m_num_labels*m_num_labels*m_num_pairs,0.0f);
}

double DiscreteCostFunction::evaluateTotalCostSumZeroLabeling() {

    int label = 0;
    double cost_sum_unary = 0.0f;
    double cost_sum_pairwise = 0.0f;
    double cost_sum_triplet = 0.0f;

    for (int i = 0; i < m_num_nodes; ++i)
        cost_sum_unary += computeUnaryCost(i, label);

    for (int p = 0; p < m_num_pairs; ++p)
        cost_sum_pairwise += computePairwiseCost(p, label, label);

    for (int t = 0; t < m_num_triplets; ++t)
        cost_sum_triplet += computeTripletCost(t, label, label, label);

    if(_debug)
        cout << "cost_sum_unary " << cost_sum_unary << " cost_sum_pairwise "
        << cost_sum_pairwise  << " cost_sum_triplet " << cost_sum_triplet
        << " total " <<  cost_sum_unary + cost_sum_pairwise + cost_sum_triplet
        << " m_num_triplets " << m_num_triplets <<   endl;

    return cost_sum_unary + cost_sum_pairwise + cost_sum_triplet;
}

double DiscreteCostFunction::evaluateTotalCostSum(const int *labeling, const int *pairs, const int *triplets, const int *quartets) {

    double cost_sum_unary = 0.0f;
    double cost_sum_pairwise = 0.0f;
    double cost_sum_triplet = 0.0f;

    for (int i = 0; i < m_num_nodes; ++i)
        cost_sum_unary += computeUnaryCost(i, labeling[i]);

    for (int p = 0; p < m_num_pairs; ++p)
        cost_sum_pairwise += computePairwiseCost(p, labeling[pairs[p * 2]], labeling[pairs[p * 2 + 1]]);

    for (int t = 0; t < m_num_triplets; ++t)
        cost_sum_triplet += computeTripletCost(t, labeling[triplets[t * 3]], labeling[triplets[t * 3 + 1]], labeling[triplets[t * 3 + 2]]);

    if(_verbosity)
        cout << "cost_sum_unary " << cost_sum_unary << " cost_sum_pairwise " << cost_sum_pairwise
        << " cost_sum_triplet " << cost_sum_triplet
        <<" total " <<  cost_sum_unary + cost_sum_pairwise + cost_sum_triplet
        << " m_num_triplets 2 " << m_num_triplets <<  endl;

    return cost_sum_unary + cost_sum_pairwise + cost_sum_triplet;
}

double DiscreteCostFunction::evaluateUnaryCostSum(const int *labeling) {

    double cost_sum_unary = 0.0f;
    for(int i = 0; i < m_num_nodes; ++i)
        cost_sum_unary += computeUnaryCost(i,labeling[i]);
    return cost_sum_unary;
}

double DiscreteCostFunction::evaluatePairwiseCostSum(const int *labeling, const int *pairs) {

    double cost_sum_pairwise = 0.0f;
    for (int p = 0; p < m_num_pairs; ++p)
        cost_sum_pairwise += computePairwiseCost(p, labeling[pairs[p * 2]], labeling[pairs[p * 2 + 1]]);
    return cost_sum_pairwise;
}

double DiscreteCostFunction::evaluateTripletCostSum(const int *labeling, const int *triplets) {

    double cost_sum_triplet = 0.0f;
    for (int t = 0; t < m_num_triplets; ++t)
        cost_sum_triplet += computeTripletCost(t, labeling[triplets[t * 3]], labeling[triplets[t * 3 + 1]],labeling[triplets[t * 3 + 2]]);
    return cost_sum_triplet;
}

//================================SURFACE CLASS===========================================================================//
void SRegDiscreteCostFunction::initialize(int numNodes,  int numLabels, int numPairs, int numTriplets)
{
    if (_TARGET.nvertices() == 0 || _SOURCE.nvertices() == 0)
        throw MeshregException("CostFunction::You must supply source and target meshes.");
    if(_HIGHREScfweight.Ncols() != _SOURCE.nvertices())
        _HIGHREScfweight.ReSize(1, _SOURCE.nvertices());
    if (_HIGHREScfweight.Nrows() != 1 && _HIGHREScfweight.Nrows() != FEAT->get_dim())
        throw MeshregException("DiscreteModel ERROR:: costfunction weighting has dimensions incompatible with data");

    DiscreteCostFunction::initialize(numNodes,numLabels,numPairs,numTriplets);
}

void SRegDiscreteCostFunction::set_parameters(myparam& ALLPARAMS){

    myparam::iterator it;
    it=ALLPARAMS.find("lambda"); _reglambda=boost::get<float>(it->second);
    it=ALLPARAMS.find("range"); _controlptrange=boost::get<float>(it->second);
    it=ALLPARAMS.find("CPres"); _RES=boost::get<int>(it->second);
    it=ALLPARAMS.find("anatres"); _aRES=boost::get<int>(it->second);
    it=ALLPARAMS.find("simmeasure"); _simmeasure=boost::get<int>(it->second); sim.set_simval(_simmeasure);
    it=ALLPARAMS.find("verbosity"); _verbosity=boost::get<bool>(it->second);
    it=ALLPARAMS.find("outdir"); _outdir=boost::get<string>(it->second);
    it=ALLPARAMS.find("regularisermode"); _rmode=boost::get<int>(it->second);
    it=ALLPARAMS.find("sigma_in"); _sigma=boost::get<float>(it->second);
    it=ALLPARAMS.find("numthreads"); _threads=boost::get<int>(it->second);
}

newresampler::Mesh SRegDiscreteCostFunction::project_anatomical() {

    newresampler::Mesh _aICOtrans = _aICO;
    barycentric_mesh_interpolation(_aICOtrans,_ORIG,_SOURCE);
    return newresampler::project_mesh(_aICOtrans,_TARGEThi,_aTARGET);
}

void SRegDiscreteCostFunction::reset_anatomical(const string &outdir, int iter) {

    NEWMAT::ColumnVector strainstmp(_aSOURCE.ntriangles());
    double perc = 0;
    _iter = iter;
    if(_aSOURCE.nvertices() > 0)
    {
        _aSOURCEtrans = project_anatomical();
        MAXstrain = 0.0;

        #pragma omp parallel for
        for (int i = 0; i < _aSOURCE.ntriangles(); i++)
            strainstmp(i + 1) = calculate_triangular_strain(i, _aSOURCE, _aSOURCEtrans, _mu, _kappa);

        for (int i = 0; i < _aSOURCE.ntriangles(); i++)
            if(strainstmp(i+1) > MAXstrain)
                MAXstrain = strainstmp(i+1);

        if(MAXstrain > 1e-8)
        {
            MISCMATHS::Histogram strainHist(strainstmp,256);
            perc = strainHist.getPercentile(0.95);
        }
        if(iter == 1) strain95 = perc;
    }
}

void SRegDiscreteCostFunction::set_initial_angles(const std::vector<std::vector<double>>& angles) {
    double meanang = 0.0;
    for (const auto& angle : angles)
        for (double j : angle)
            meanang += j;

    _MEANANGLE = meanang / (3.0*angles.size());
}

bool SRegDiscreteCostFunction::within_controlpt_range(int CPindex, int sourceindex) {

    newresampler::Point CP = _CPgrid.get_coord(CPindex);
    newresampler::Point SP = _SOURCE.get_coord(sourceindex);
    double dist = 2 * RAD * asin((CP-SP).norm()/(2*RAD));
    return (dist < _controlptrange * MAXSEP(CPindex + 1));
}
//================================Non Linear SURFACE CLASS===========================================================================//
NonLinearSRegDiscreteCostFunction::NonLinearSRegDiscreteCostFunction() {
    _expscaling = 1; _k_exp = 2.0; _maxdist = 4; _rexp = 2; _kNN = 5; _rmode = 1; _mu = 0.4; _kappa = 1.6;
}

void NonLinearSRegDiscreteCostFunction::initialize(int numNodes, int numLabels, int numPairs, int numTriplets) {
    SRegDiscreteCostFunction::initialize(numNodes, numLabels, numPairs, numTriplets);
}

void NonLinearSRegDiscreteCostFunction::set_parameters(myparam& ALLPARAMS) {

    myparam::iterator it;
    it=ALLPARAMS.find("maxdistortion");_maxdist=boost::get<float>(it->second);
    it=ALLPARAMS.find("exponent");_rexp=boost::get<float>(it->second);
    it=ALLPARAMS.find("weight");_dweight=boost::get<bool>(it->second);
    it=ALLPARAMS.find("anorm");_anorm=boost::get<bool>(it->second);
    it=ALLPARAMS.find("scaling");_expscaling=boost::get<float>(it->second);
    it=ALLPARAMS.find("shearmodulus");_mu=boost::get<float>(it->second);
    it=ALLPARAMS.find("bulkmodulus");_kappa=boost::get<float>(it->second);
    it=ALLPARAMS.find("kexponent");_k_exp=boost::get<float>(it->second);
    it=ALLPARAMS.find("kNN");_kNN=boost::get<int>(it->second);
    it=ALLPARAMS.find("pottsthreshold");_pottsthreshold=boost::get<float>(it->second);

    SRegDiscreteCostFunction::set_parameters(ALLPARAMS);
}

double NonLinearSRegDiscreteCostFunction::computeTripletCost(int triplet, int labelA, int labelB, int labelC){

    if (!_triplets) throw MeshregException("DiscreteModel ERROR:: must run settripletss() prior to computTripleCost ");

    if(triplet == 0 && _debug) { sumlikelihood = 0; sumregcost = 0; }

    double cost = 0.0, weight = 1.0;

    std::vector<int> id;
    std::map<int,newresampler::Point> vertex;

    id.push_back(_triplets[3*triplet]);
    id.push_back(_triplets[3*triplet+1]);
    id.push_back(_triplets[3*triplet+2]);

    vertex[id[0]] = (*ROTATIONS)[id[0]] * _labels[labelA];
    vertex[id[1]] = (*ROTATIONS)[id[1]] * _labels[labelB];
    vertex[id[2]] = (*ROTATIONS)[id[2]] * _labels[labelC];

    double likelihood = triplet_likelihood(triplet, id[0], id[1], id[2], vertex[id[0]], vertex[id[1]], vertex[id[2]]);

    newresampler::Triangle TRI(vertex[id[0]], vertex[id[1]], vertex[id[2]], 0);
    newresampler::Triangle TRI_noDEF(_CPgrid.get_coord(_triplets[3*triplet]),
                                     _CPgrid.get_coord(_triplets[3*triplet+1]),
                                     _CPgrid.get_coord(_triplets[3*triplet+2]),0);

    // only estimate cost if it doesn't cause folding
    if((TRI.normal() | TRI_noDEF.normal()) < 0)
    {
        cost = 1.0;
        weight += 1e6;
    }
    else
    {
        switch(_rmode)
        {
            case 2:
            {
                double tweight = 1.0, diff = 0.0;
                newresampler::Triangle TRI_ORIG(_ORIG.get_coord(_triplets[3*triplet]),
                                                _ORIG.get_coord(_triplets[3*triplet+1]),
                                                _ORIG.get_coord(_triplets[3*triplet+2]), 0);
                vector<double> deformed_angles = TRI.get_angles();
                vector<double> orig_angles = TRI_ORIG.get_angles();
                double distortion = log2(TRI.get_area()/TRI_ORIG.get_area());

                if (_dweight && distortion != 1) tweight = exp(abs(distortion - 1));
                else tweight = 1.0;

                for (int i = 0; i < 3; i++)
                    diff += std::pow(deformed_angles[i] - orig_angles[i], 2);

                cost = tweight*(sqrt(diff));

                if (_anorm) weight *= (1 / _MEANANGLE);
                break;
            }
            case 3:
            {
                newresampler::Triangle TRI_ORIG(_ORIG.get_coord(_triplets[3*triplet]),
                                                _ORIG.get_coord(_triplets[3*triplet+1]),
                                                _ORIG.get_coord(_triplets[3*triplet+2]), 0);

                cost = calculate_triangular_strain(TRI_ORIG, TRI, _mu, _kappa, std::shared_ptr<NEWMAT::ColumnVector>(), _k_exp);
                break;
            }
            case 4:
            {
                double diff = 0.0;
                map<int,bool> moved2;
                map<int,newresampler::Point> transformed_points;

                for (unsigned int n = 0; n < NEARESTFACES[triplet].size(); n++)
                {
                    newresampler::Triangle TRIorig = _aSOURCE.get_triangle(NEARESTFACES[triplet][n]);
                    newresampler::Triangle TRItrans = deform_anatomy(triplet, n, vertex, moved2, transformed_points);
                    vector<double> deformed_angles = TRItrans.get_angles();
                    vector<double> orig_angles = TRIorig.get_angles();
                    double distortion = log2(TRItrans.get_area()/TRIorig.get_area());
                    double tweight = 0.0;
                    if (_dweight && distortion != 1)
                        tweight = exp(abs(distortion - 1));
                    else
                        tweight = 1;
                    for (int i = 0; i < 3; i++)
                        diff += std::pow(deformed_angles[i] - orig_angles[i], 2);
                    diff = tweight * sqrt(diff);
                    cost += diff;
                }
                cost = cost / (double) NEARESTFACES[triplet].size();
                break;
            }
            case 5:
            {
                map<int,bool> moved2;
                map<int,newresampler::Point> transformed_points;
                for (unsigned int n = 0; n < NEARESTFACES[triplet].size(); n++)
                {
                    newresampler::Triangle TRIorig = _aSOURCE.get_triangle(NEARESTFACES[triplet][n]);
                    newresampler::Triangle TRItrans = deform_anatomy(triplet, n, vertex, moved2, transformed_points);
                    cost += calculate_triangular_strain(TRIorig, TRItrans, _mu, _kappa, std::shared_ptr<NEWMAT::ColumnVector>(), _k_exp);
                }
                cost = cost / NEARESTFACES[triplet].size();
                break;
            }
            default:
                throw MeshregException("DiscreteModel computeTripletCost regoption does not exist");
        }
    }

    if (abs(cost) < 1e-8) cost = 0.0;

    if(_debug)
    {
        sumlikelihood += likelihood;
        sumregcost += weight * _reglambda * MISCMATHS::pow(cost,_rexp);
    }
    cost = likelihood + weight * _reglambda * MISCMATHS::pow(cost,_rexp);

    return cost; // normalise to try and ensure equivalent lambda for each resolution level
}

double NonLinearSRegDiscreteCostFunction::computePairwiseCost(int pair, int labelA, int labelB){

    if (!_pairs) throw MeshregException("DiscreteModel ERROR:: must run setPairs() prior to computePairwiseCost ");

    NEWMAT::Matrix R1 = estimate_rotation_matrix(_CPgrid.get_coord(_pairs[2*pair]),(*ROTATIONS)[_pairs[2*pair]]*_labels[labelA]);
    NEWMAT::Matrix R2 = estimate_rotation_matrix(_CPgrid.get_coord(_pairs[2*pair+1]),(*ROTATIONS)[_pairs[2*pair+1]]*_labels[labelB]);
    NEWMAT::Matrix R_diff = R1.t() * R2;

    const double theta_MVD = 2 * asin(MVDmax/(2*RAD));
    const double theta = acos((R_diff.Trace()-1)/2);
    double weight = 1.0, totalcost = 0.0, cost = 0.0;

    newresampler::Point v0 = _CPgrid.get_coord(_pairs[2*pair]);
    newresampler::Point v1 = _CPgrid.get_coord(_pairs[2*pair+1]);

    _CPgrid.set_coord(_pairs[2*pair],(*ROTATIONS)[_pairs[2*pair]]*_labels[labelA]);
    _CPgrid.set_coord(_pairs[2*pair+1],(*ROTATIONS)[_pairs[2*pair+1]]*_labels[labelB]);

    if(fabs(1 - (R_diff.Trace() - 1) / 2) > EPSILON)
    {
        for(auto j = _CPgrid.tIDbegin(_pairs[2*pair]); j != _CPgrid.tIDend(_pairs[2*pair]); j++)
        {
            if((_oCPgrid.get_triangle(*j).normal() | _CPgrid.get_triangle(*j).normal()) < 0)
            {
                cost = 1.0;
                weight += 1e6;
                break;
            }
        }

        cost = (sqrt(2) * theta) / theta_MVD;

        if(_rexp == 1)
            totalcost = weight * _reglambda * cost;
        else
            totalcost = weight * _reglambda * MISCMATHS::pow(cost,_rexp);
    }

    _CPgrid.set_coord(_pairs[2*pair],v0);
    _CPgrid.set_coord(_pairs[2*pair+1],v1);

    return totalcost;
}

void NonLinearSRegDiscreteCostFunction::computePairwiseCosts(const int *pairs) {

    for (int i = 0; i < m_num_pairs; i++)
        for (unsigned int j = 0; j < _labels.size(); j++)
            for (unsigned int k = 0; k < _labels.size(); k++)
                paircosts[i * m_num_labels * m_num_labels + k * m_num_labels + j] = computePairwiseCost(i, j, k);
}

void  NonLinearSRegDiscreteCostFunction::computeUnaryCosts(){
    // for each control point resample data into blocks each influencing a single control point
    // calculates similarity
    for (int j = 0; j < m_num_labels; j++)
        for (int k = 0; k < _CPgrid.nvertices(); k++)
            unarycosts[j * m_num_nodes + k] = computeUnaryCost(k, j);
}

newresampler::Triangle NonLinearSRegDiscreteCostFunction::deform_anatomy(int trip, int n, map<int,newresampler::Point>& vertex,
                                                                    map<int,bool>& moved, map<int,newresampler::Point>& transformed) {
    newresampler::Triangle TRItrans;
    newresampler::Point newPt;
    int tindex;
    map<int,double> weight;

    for(int i=0;i<3;i++)
    { // for each point in face
        tindex = _aSOURCE.get_triangle(NEARESTFACES[trip][n]).get_vertex_no(i);

        if(moved.find(tindex) == moved.end())
        {
            moved[tindex] = true;

            newPt.X = 0; newPt.Y = 0; newPt.Z = 0;
            for (auto &it: _ANATbaryweights[tindex])
                newPt += vertex[it.first] * it.second;

            newresampler::Triangle closest_triangle = anattree->get_closest_triangle(newPt);

            if(closest_triangle.get_no() == -1) {
                std::cout << newPt << std::endl;
                throw MeshregException("closest triangle not found");
            }

            newresampler::Point v0 = closest_triangle.get_vertex_coord(0),
                                v1 = closest_triangle.get_vertex_coord(1),
                                v2 = closest_triangle.get_vertex_coord(2);
            int n0 = closest_triangle.get_vertex_no(0),
                n1 = closest_triangle.get_vertex_no(1),
                n2 = closest_triangle.get_vertex_no(2);

            weight = newresampler::calc_barycentric_weights(v0,v1,v2,newPt,n0,n1,n2);

            newPt.X = 0; newPt.Y = 0; newPt.Z = 0;
            for (auto &it: weight)
                newPt += _aTARGET.get_coord(it.first) * it.second;

            transformed[tindex]=newPt;
            TRItrans.set_vertex(i, newPt);
        }
        else
            TRItrans.set_vertex(i, transformed[tindex]);
    }

    return TRItrans;
}

void NonLinearSRegDiscreteCostFunction::resample_weights(){
// TAKE DATA TERM WEIGHTING AS MAXIMUM OF WEIGHTS ACROSS ALL DIMENSIONS

    AbsoluteWeights.ReSize(_SOURCE.nvertices());
    AbsoluteWeights = 0;

    #pragma omp parallel for
    for (int k = 1; k <= _SOURCE.nvertices(); k++)
    {
        double maxweight = std::numeric_limits<double>::lowest();
        for (int j = 1; j <= _HIGHREScfweight.Nrows(); j++)
            if (_HIGHREScfweight(j, k) > maxweight)
                maxweight = _HIGHREScfweight(j, k);

        AbsoluteWeights(k) = maxweight;
    }

    newresampler::Mesh tmp = _SOURCE;
    tmp.set_pvalues(AbsoluteWeights);
    AbsoluteWeights = newresampler::metric_resample(tmp, _CPgrid).get_pvalues();
}

void NonLinearSRegDiscreteCostFunction::get_target_data(int node, const NEWMAT::Matrix& PtROTATOR) {

    _targetdata[node].clear();
    _targetdata[node].resize(_sourceinrange[node].size());

    #pragma omp parallel for
    for(unsigned int i = 0; i < _sourceinrange[node].size(); i++)
    {
        newresampler::Point tmp = PtROTATOR * _SOURCE.get_coord(_sourceinrange[node][i]);

        newresampler::Triangle closest_triangle = targettree->get_closest_triangle(tmp);

        newresampler::Point v0 = closest_triangle.get_vertex_coord(0),
                            v1 = closest_triangle.get_vertex_coord(1),
                            v2 = closest_triangle.get_vertex_coord(2);
        int n0 = closest_triangle.get_vertex_no(0),
            n1 = closest_triangle.get_vertex_no(1),
            n2 = closest_triangle.get_vertex_no(2);

        _targetdata[node][i] = newresampler::barycentric_weight(v0, v1, v2, tmp,
                                                                         FEAT->get_ref_val(1, n0 + 1),
                                                                         FEAT->get_ref_val(1, n1 + 1),
                                                                         FEAT->get_ref_val(1, n2 + 1));
    }
}

//================================ UNIVARIATE Non Linear SURFACE CLASS===========================================================================//
void UnivariateNonLinearSRegDiscreteCostFunction::initialize(int numNodes, int numLabels, int numPairs, int numTriplets) {

    NonLinearSRegDiscreteCostFunction::initialize(numNodes, numLabels, numPairs, numTriplets);
    _sourcedata.clear(); _sourcedata.resize(_CPgrid.nvertices(), vector<double>());
    _sourceinrange.clear(); _sourceinrange.resize(_CPgrid.nvertices(), vector<int>());
    _targetdata.clear(); _targetdata.resize(_CPgrid.nvertices(), vector<double>());
    _weights.clear(); _weights.resize(_CPgrid.nvertices(), vector<double>());
}

void UnivariateNonLinearSRegDiscreteCostFunction::get_source_data() {

    for (auto& i : _sourcedata) i.clear();

    #pragma omp parallel for
    for (int k = 0; k < _CPgrid.nvertices(); k++)
        for (int i = 0; i < _SOURCE.nvertices(); i++)
            if (within_controlpt_range(k, i))
            {
                _sourceinrange[k].push_back(i);
                _sourcedata[k].push_back(FEAT->get_input_val(1, i + 1));
                if (_HIGHREScfweight.Nrows() >= 1)
                    _weights[k].emplace_back(_HIGHREScfweight(1, i + 1));
                else
                    _weights[k].emplace_back(1.0);
            }
    resample_weights();
}

double UnivariateNonLinearSRegDiscreteCostFunction::computeUnaryCost(int node, int label){

    NEWMAT::Matrix ROT = estimate_rotation_matrix(_CPgrid.get_coord(node),(*ROTATIONS)[node] * _labels[label]); // this Matrix will rotate control point k to the label

    get_target_data(node,ROT);

    double cost = AbsoluteWeights(node + 1) *
            sim.get_sim_for_min(_sourcedata[node], _targetdata[node],_weights[node]);// weights similarity by cost function weighting

    if(_simmeasure == 1 || _simmeasure == 2)
        return cost;
    else
        return -cost;
}

//================================Multivariate Non Linear SURFACE CLASS===========================================================================//
void MultivariateNonLinearSRegDiscreteCostFunction::initialize(int numNodes, int numLabels, int numPairs, int numTriplets)
{
    NonLinearSRegDiscreteCostFunction::initialize(numNodes, numLabels,numPairs,numTriplets);
    _sourcedata.clear(); _sourcedata.resize(_CPgrid.nvertices(), vector<double> ());
    _sourceinrange.clear(); _sourceinrange.resize(_CPgrid.nvertices(), vector<int> ());
    _targetdata.clear(); _targetdata.resize(_CPgrid.nvertices(), vector<double> ());
    _weights.clear(); _weights.resize(_CPgrid.nvertices(), vector<double> ());
}

void MultivariateNonLinearSRegDiscreteCostFunction::get_source_data() {

    for(auto& i : _sourcedata) i.clear();

    #pragma omp parallel for
    for (int k = 0; k < _CPgrid.nvertices(); k++)
        for (int i = 0; i < _SOURCE.nvertices(); i++)
            if (within_controlpt_range(k, i))
            {
                _sourceinrange[k].push_back(i);
                for (int d = 1; d <= FEAT->get_dim(); d++)
                {
                    _sourcedata[k].push_back(FEAT->get_input_val(d, k + 1));
                    if (_HIGHREScfweight.Nrows() >= d)
                        _weights[k].emplace_back(_HIGHREScfweight(d, k + 1));
                    else
                        _weights[k].emplace_back(1.0);
                }
            }
    resample_weights();
}

double MultivariateNonLinearSRegDiscreteCostFunction::computeUnaryCost(int node, int label){

    double cost	= 0.0;

    get_target_data(node,estimate_rotation_matrix(_CPgrid.get_coord(node),(*ROTATIONS)[node] * _labels[label]));

    for(int i = 0; i < _sourceinrange[node].size(); ++i)
        for(int dim = 0; dim < FEAT->get_dim(); ++dim)
        {
            int offset = dim * FEAT->get_dim();
            cost += sim.get_sim_for_min({_sourcedata[node].begin() + offset, _sourcedata[node].begin() + offset+FEAT->get_dim()-1 },
                                        {_targetdata[node].begin() + offset, _targetdata[node].begin() + offset+FEAT->get_dim()-1 },
                                        {_weights[node].begin() + offset,_weights[node].begin() + offset+FEAT->get_dim()-1 });
        }

    if(!_sourceinrange[node].empty()) cost /= (_sourceinrange[node].size() * FEAT->get_dim());

    if(_simmeasure==1 || _simmeasure==2)
        return AbsoluteWeights(node + 1) * cost;
    else
        return -AbsoluteWeights(node + 1) * cost;
}

} //namespace newmeshreg
