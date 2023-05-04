#ifndef NEWMESHREG_DISCRETEMODEL_H
#define NEWMESHREG_DISCRETEMODEL_H

#include <omp.h>

#include "DiscreteCostFunction.h"

namespace newmeshreg {

class DiscreteModel {

public:
    DiscreteModel()
        : m_num_nodes(0), m_num_labels(0), m_num_pairs(0), m_num_triplets(0),
        labeling(nullptr), pairs(nullptr), triplets(nullptr), m_verbosity(false) {}

    explicit DiscreteModel(myparam& P)
        : m_num_nodes(0), m_num_labels(0), m_num_pairs(0), m_num_triplets(0),
        labeling(nullptr), pairs(nullptr), triplets(nullptr), m_verbosity(false) {}

    virtual ~DiscreteModel() {
        delete[] labeling; labeling = nullptr;
        delete[] pairs; pairs = nullptr;
        delete[] triplets; triplets = nullptr;
    }

    //---GET---//
    int getNumNodes() const { return m_num_nodes; }
    int getNumLabels() const { return m_num_labels; }
    int getNumPairs() const { return m_num_pairs; }
    int getNumTriplets() const { return m_num_triplets; }
    int* getLabeling() { return labeling; }
    const int* getLabeling() const { return labeling; }
    const int* getPairs() const { return pairs; }
    const int* getTriplets() const { return triplets; }
    virtual std::shared_ptr<DiscreteCostFunction> getCostFunction() = 0;

    //---COMPUTE---//
    virtual void computeUnaryCosts(){}
    virtual double computeUnaryCost(int node, int label) { return 0; }
    virtual void computePairwiseCosts() {}
    virtual double computePairwiseCost(int pair, int labelA, int labelB) { return 0; }
    virtual double computeTripletCost(int triplet, int labelA, int labelB, int labelC) { return 0; }
    virtual double evaluateTotalCostSum() { return 0; }

    //---MODIFY---//
    virtual void applyLabeling(){}
    virtual void applyLabeling(int *discreteLabeling) {}
    virtual void report(){}

    //---INITIALIZE---//
    void resetLabeling() { if(labeling) std::fill(labeling,labeling + m_num_nodes,0.0); }
    virtual void Initialize(const newresampler::Mesh&) {}
    virtual void setupCostFunction() {}
    virtual void set_parameters(myparam& PAR) {}

protected:
    void initLabeling() {
        if(m_num_nodes != 0)
        {
            delete[] labeling;
            labeling = new int[m_num_nodes];
            resetLabeling();
        }
    }

    int m_num_nodes;    // Number of model nodes (e.g. grid nodes, variables, etc).
    int m_num_labels;   // Number of labels.
    int m_num_pairs;    // Number of node pairs.
    int m_num_triplets; // Number of node triplets.
    int* labeling;      // Labeling array.
    int* pairs;         // Node pairs array.
    int* triplets;      // Node triplets array.

    std::string m_outdir;
    bool m_verbosity;
    int _nthreads = 1;
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

class SRegDiscreteModel : public DiscreteModel {

public:
    SRegDiscreteModel() = default;

    explicit SRegDiscreteModel(myparam& PAR) {
        set_parameters(PAR);
        initialize_cost_function(m_multivariate, PAR);
    }

    ~SRegDiscreteModel() override = default;

    std::shared_ptr<DiscreteCostFunction> getCostFunction() override { // upcast
        std::shared_ptr<DiscreteCostFunction> dcostfct = costfct;
        return dcostfct;
    }

    void computeUnaryCosts() override { costfct->computeUnaryCosts(); }
    double computeUnaryCost(int node, int label) override { return (costfct) ? costfct->computeUnaryCost(node,label) : 0.0f; }
    void computePairwiseCosts() override { costfct->computePairwiseCosts(pairs); }
    double computePairwiseCost(int pair, int labelA, int labelB) override {return (costfct) ? costfct->computePairwiseCost(pair,labelA,labelB): 0.0f; }
    double computeTripletCost(int triplet, int labelA, int labelB, int labelC) override { return (costfct) ? costfct->computeTripletCost(triplet,labelA,labelB,labelC) : 0.0f; }
    double evaluateTotalCostSum() override { return (costfct) ? costfct->evaluateTotalCostSum(labeling,pairs,triplets) : 0.0f; }
    void set_parameters(myparam& PAR) override;

    //---INITIALIZE MODEL---//
    virtual void set_meshspace(const newresampler::Mesh& target, const newresampler::Mesh& source) { m_TARGET = target; m_SOURCE = source; }
    void set_anatomical_meshspace(const newresampler::Mesh& ref_sphere, const newresampler::Mesh& ref_anat,
                                  const newresampler::Mesh& source_sphere, const newresampler::Mesh & source_anat) {
        costfct->set_anatomical(ref_sphere, ref_anat, source_sphere, source_anat);
    }
    void set_anatomical_neighbourhood(const std::vector<std::map<int,double>>& weights, const std::vector<std::vector<int>>& neighbourhood) {
        costfct->set_anatomical_neighbourhood(weights, neighbourhood);
    }

    void set_featurespace(const std::shared_ptr<featurespace>& FEATURES) {
        costfct->set_featurespace(FEATURES);
    }

    // costfunction weighting combines source and reference weightings at beginning of optimisation iteration -
    //will not be 100% accurate but will remove any sensitivity of the label choices to weighting
    void setupCostFunctionWeighting(const NEWMAT::Matrix& Weight) { costfct->set_dataaffintyweighting(Weight); }

    // source needs to be reset after every iteration of discrete optimisation
    virtual void reset_meshspace(const newresampler::Mesh& source) {
        m_SOURCE = source;
        costfct->reset_source(source);
    }

    virtual void reset_CPgrid(const newresampler::Mesh& grid) { m_CPgrid = grid; }
    virtual void warp_CPgrid(newresampler::Mesh& START, newresampler::Mesh& END) {
        barycentric_mesh_interpolation(m_CPgrid,START,END, _nthreads);
        unfold(m_CPgrid);
    }

    void initialize_cost_function(bool MV, myparam &P);
    void Initialize_sampling_grid();
    void label_sampling_grid(int, double, newresampler::Mesh&);
    std::vector<newresampler::Point> rescale_sampling_grid();

    void Initialize(const newresampler::Mesh&) override;

    virtual void set_debug(){ m_debug = true; costfct->debug(); } // for debuging

    newresampler::Mesh get_TARGET(){ return m_TARGET; }
    virtual newresampler::Mesh get_CPgrid() { return m_CPgrid; }

    void report() override { if (costfct) costfct->report(); }

protected:
    newresampler::Mesh m_TARGET; // TARGET MESH
    newresampler::Mesh m_SOURCE; // SOURCE MESH
    newresampler::Mesh m_CPgrid; // CONTROL POINT GRID
    newresampler::Mesh m_samplinggrid;

    std::shared_ptr<newresampler::Octree> m_inputtree;

    int m_SGres = 4; // sampling grid resolution
    int m_iter = 0; // iteration of the discrete optimisation
    int m_centroid = 0; // used for selecting which sampling grid vertex will form the center of the sampling grid
    int m_regoption = 2; // sim measure i.e. correlation
    double m_maxs_dist = 0.0; // define maximum distance between the centre of the sampling grid and the furthest label
    double MVD = 0.0;
    float m_scale = 0.0;
    bool m_multivariate = false;
    bool m_debug = false;
    bool m_triclique = false;
    bool _pairwise = false;
    bool m_rescalelabels = false;
    newresampler::Point centre;

    std::vector<newresampler::Point> m_samples;  // samples based on  vertices of sampling grid
    std::vector<newresampler::Point> m_barycentres; // samples based on barycentres of sampling grid
    std::vector<newresampler::Point> m_labels; // labels iterates between samples and barycnetres and is the label set used within cosfct
    std::vector<NEWMAT::Matrix> m_ROT; // rotates sampling grid to each control point
    std::shared_ptr<SRegDiscreteCostFunction> costfct;  // costfunction object
};

class NonLinearSRegDiscreteModel: public SRegDiscreteModel {
public:
    explicit NonLinearSRegDiscreteModel(myparam& P) : SRegDiscreteModel(P){  set_parameters(P); }

    void applyLabeling() override { applyLabeling(labeling); }
    void applyLabeling(int* discreteLabeling) override;

    void estimate_pairs();
    void estimate_triplets();

    void Initialize(const newresampler::Mesh&) override;

    void get_rotations(std::vector<NEWMAT::Matrix>&);
    void setupCostFunction() override;
};

} //namespace newmeshreg

#endif //NEWMESHREG_DISCRETEMODEL_H
