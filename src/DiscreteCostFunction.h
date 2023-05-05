#ifndef NEWMESHREG_DISCRETECOSTFUNCTION_H
#define NEWMESHREG_DISCRETECOSTFUNCTION_H

#include <omp.h>
#include "featurespace.h"
#include "newresampler/resampler.h"
#include "similarities.h"

#define RAD 100.0
#define EPSILON 1.0E-8

namespace newmeshreg {

class DiscreteCostFunction {

public:
    DiscreteCostFunction() = default;

    virtual ~DiscreteCostFunction(){
        delete[] unarycosts;
        delete[] paircosts;
    }

    //---SET--//
    virtual void reset(); //Resets all costs to zero.
    void setPairs(int* p) { _pairs = p; } //Sets the pairs index buffer.
    void setTriplets(int* p) { _triplets = p; } //Sets the pairs index buffer.

    //---GET--//
    double* getUnaryCosts() { return unarycosts; } //Returns the unary costs look-up table.
    double* getPairwiseCosts() { return paircosts; } //Returns the pairwise costs look-up table.

    //---COMPUTE--//
    virtual void computeUnaryCosts() {}; //Computes the unary costs look-up table.
    virtual double computeUnaryCost(int node, int label) { return 0; } //Computes the unary potential for the given node

    virtual void computePairwiseCosts(const int *pairs) {} //Computes the pairwise costs look-up table.
    virtual double computePairwiseCost(int pair, int labelA, int labelB) { return 0; } //Computes the pairwise potential for a the given pair and labels.

    virtual double computeTripletCost(int triplet, int labelA, int labelB, int labelC) { return 0; } //Computes the triplet potential for a the given triplet and labels.
    virtual double evaluateTotalCostSum(const int *labeling, const int *pairs, const int *triplets); //Evaluates the total cost for the given labeling.

    virtual void set_parameters(myparam&) = 0;
    virtual void report() {}

protected:
    virtual void initialize(int numNodes, int numLabels, int numPairs, int numTriplets);

    int m_num_nodes = 0;
    int	m_num_labels = 0;
    int m_num_pairs = 0;
    int m_num_triplets = 0;

    double* unarycosts = nullptr; // Unary potentials look-up table.
    double* paircosts = nullptr; // Pairwise potentials look-up table.

    int* _pairs = nullptr;
    int* _triplets = nullptr;

    float _reglambda = 0.0;  // scaling parameter for regulariser

    bool _debug = false;
    bool _verbosity = false;
};

// should be able to implement affine spherical registration using discrete labels represent fixed rotations
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

    void set_parameters(myparam &) override {}

protected:
    std::map<int,std::vector<double>> unaryenergies; // maps of nodes  xlabels x vals
    std::map<int,std::vector<double>> pairenergies;
};

class SRegDiscreteCostFunction: public DiscreteCostFunction {

public:
    //---INITIALISE---//
    SRegDiscreteCostFunction() = default;

    void set_parameters(myparam&) override;
    void initialize(int numNodes, int numLabels, int numPairs, int numTriplets) override;

    void initialize_regulariser() {
        if (_aSOURCE.nvertices() > 0 && _rmode >= 3)
            anattree = std::make_shared<newresampler::Octree>(_TARGEThi);
    }

    //---SET---//
    // data input and reference mesh, plus low resolution control point grid
    virtual void set_meshes(const newresampler::Mesh& target,const newresampler::Mesh& source, const newresampler::Mesh& GRID){
        _TARGET = target; _SOURCE = source; _ORIG = source; _CPgrid = GRID; _oCPgrid = GRID;
    }
    void set_anatomical(const newresampler::Mesh& targetS, const newresampler::Mesh& targetA, const newresampler::Mesh& sourceS, const newresampler::Mesh& sourceA) {
        _TARGEThi = targetS; _aTARGET = targetA; _aICO = sourceS, _aSOURCE = sourceA;
    }
    void set_anatomical_neighbourhood(const std::vector<std::map<int,double>>& weights, const std::vector<std::vector<int>>& neighbourhood) {
        _ANATbaryweights = weights; NEARESTFACES = neighbourhood;
    }
    void set_featurespace(const std::shared_ptr<featurespace>& features) {
        FEAT = features;
    }
    void set_labels(const std::vector<newresampler::Point>& labellist, const std::vector<NEWMAT::Matrix>& ROT = std::vector<NEWMAT::Matrix>()) {
        _labels = labellist;
        ROTATIONS = std::make_shared<std::vector<NEWMAT::Matrix>>(ROT);
    }
    virtual void set_spacings(const NEWMAT::ColumnVector& spacings, double MAX) { MAXSEP = spacings; MVDmax = MAX; }
    void set_dataaffintyweighting(const NEWMAT::Matrix& HRWeight) { _HIGHREScfweight = HRWeight; }
    void set_octrees(std::shared_ptr<newresampler::Octree>& targett) { targettree = targett; }

    virtual void reset_source(const newresampler::Mesh& source) { _SOURCE = source; }
    virtual void reset_CPgrid(const newresampler::Mesh& grid) { _CPgrid = grid; }
    void reset_anatomical();
    void set_initial_angles(const std::vector<std::vector<double>>& angles);

    void report() override { if(_debug) std::cout << " sumlikelihood " << sumlikelihood << " sumregcost " << sumregcost <<std::endl; }
    void debug() { _debug = true; } // for debuging

    newresampler::Mesh project_anatomical();
    bool within_controlpt_range(int CPindex, int sourceindex);

    virtual void resample_weights() {}
    virtual void get_source_data() {}
    virtual double triplet_likelihood(int, int, int, int, const newresampler::Point&, const newresampler::Point&, const newresampler::Point &) { return 0; }

protected:
    //---MESHES---//
    newresampler::Mesh _TARGET; // TARGET MESH
    newresampler::Mesh _TARGEThi; // TARGET MESH
    newresampler::Mesh _SOURCE; // SOURCE MESH
    newresampler::Mesh _aTARGET; // ANATOMICAL TARGET MESH
    newresampler::Mesh _aSOURCE; // ANATOMICAL SOURCE MESH
    newresampler::Mesh _aICO; // ANATOMICAL SOURCE MESH
    newresampler::Mesh _aSOURCEtrans; // ANATOMICAL SOURCE MESH
    newresampler::Mesh _ORIG; // NON DEFORMED SOURCE
    newresampler::Mesh _oCPgrid; // NON DEFORMED CP GRID
    newresampler::Mesh _CPgrid; // CONTRO

    //---CF WEIGHTINGS---//
    NEWMAT::Matrix _HIGHREScfweight; // source mesh cfweight
    NEWMAT::ColumnVector MAXSEP;

    //---NEIGHBOURHOOD INFO---//
    std::shared_ptr<newresampler::Octree> targettree;
    std::shared_ptr<newresampler::Octree> anattree;

    std::vector<std::vector<int>> _sourceinrange;
    std::vector<std::vector<int>> NEARESTFACES;

    //---FEATURE SET---//
    std::shared_ptr<featurespace> FEAT; // holds data
    sparsesimkernel sim; // similarity object

    //---LABELLING PARAMETERS---//
    std::vector<newresampler::Point> _labels;
    std::shared_ptr<std::vector<NEWMAT::Matrix>> ROTATIONS; // rotates label set onto each control point

    double MVDmax = 0.0; // max distance between CPs
    double _MEANANGLE = 0.0;
    float _controlptrange = 1.0;

    double sumlikelihood = 0.0;
    double sumregcost = 0.0;

    float _mu = 0.0; // shear modulus
    float _kappa = 0.0; // bulk modulus
    float _rexp = 0.0;

    //---USER DEFINED PARAMETERS---//
    int _simmeasure = 2;
    int  _rmode = 1;
    int _threads = 1;
    float _k_exp = 0.0;
    double MAXstrain = 0.0;

    std::vector<std::vector<double>> _sourcedata;
    std::vector<std::vector<double>> _targetdata;
    std::vector<std::vector<double>> _weights;
    std::vector<std::map<int,double>> _ANATbaryweights;
};

class NonLinearSRegDiscreteCostFunction: public SRegDiscreteCostFunction {
protected:
    //---REGULARISER OPTIONS---//
    bool _dweight = false;
    bool _anorm = false;

    NEWMAT::ColumnVector AbsoluteWeights;

public:
    NonLinearSRegDiscreteCostFunction();

    void initialize(int numNodes, int numLabels, int numPairs, int numTriplets) override;
    void set_parameters(myparam&) override;

    void computeUnaryCosts() override;

    double computePairwiseCost(int pair, int labelA, int labelB) override;
    void computePairwiseCosts(const int *pairs) override;

    double computeTripletCost(int triplet, int labelA, int labelB, int labelC) override;
    double triplet_likelihood(int, int, int, int, const newresampler::Point&, const newresampler::Point&, const newresampler::Point&) override { return 0; }

    newresampler::Triangle deform_anatomy(int, int, std::map<int,newresampler::Point>&, std::map<int,bool>&, std::map<int,newresampler::Point>&);
    void resample_weights() override;

    virtual void get_target_data(int, const NEWMAT::Matrix&) = 0;
};

class UnivariateNonLinearSRegDiscreteCostFunction: public NonLinearSRegDiscreteCostFunction {
public:
    UnivariateNonLinearSRegDiscreteCostFunction() = default;
    void initialize(int numNodes, int numLabels, int numPairs, int numTriplets) override;
    void get_source_data() override;
    void get_target_data(int node, const NEWMAT::Matrix& PtROTATOR) override;
    double computeUnaryCost(int node, int label) override;
};

class MultivariateNonLinearSRegDiscreteCostFunction: public NonLinearSRegDiscreteCostFunction {
public:
    MultivariateNonLinearSRegDiscreteCostFunction() = default;
    void initialize(int numNodes, int numLabels, int numPairs, int numTriplets) override;
    void get_source_data() override;
    void get_target_data(int node, const NEWMAT::Matrix& PtROTATOR) override;
    double computeUnaryCost(int node, int label) override;
};

class HOUnivariateNonLinearSRegDiscreteCostFunction: public UnivariateNonLinearSRegDiscreteCostFunction {
public:
    HOUnivariateNonLinearSRegDiscreteCostFunction() = default;
    void initialize(int numNodes, int numLabels, int numPairs, int numTriplets) override;
    void get_source_data() override;
    void get_target_data(int triplet, const newresampler::Point& new_CP0, const newresampler::Point& new_CP1, const newresampler::Point& new_CP2);
    inline double computeUnaryCost(int node, int label) override { return 0; }
    double triplet_likelihood(int, int, int, int, const newresampler::Point&, const newresampler::Point&, const newresampler::Point&) override;
};

class HOMultivariateNonLinearSRegDiscreteCostFunction: public MultivariateNonLinearSRegDiscreteCostFunction {
public:
    HOMultivariateNonLinearSRegDiscreteCostFunction() = default;
    void initialize(int numNodes, int numLabels, int numPairs, int numTriplets) override;
    void get_source_data() override;
    void get_target_data(int triplet, const newresampler::Point& new_CP0, const newresampler::Point& new_CP1, const newresampler::Point& new_CP2);
    inline double computeUnaryCost(int node, int label) override { return 0; }
    double triplet_likelihood(int, int, int, int, const newresampler::Point&, const newresampler::Point&, const newresampler::Point&) override;
};

} //namespace newmeshreg

#endif //NEWMESHREG_DISCRETECOSTFUNCTION_H
