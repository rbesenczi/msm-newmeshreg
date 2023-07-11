#ifndef NEWMESHREG_DISCRETEGROUPCOSTFUNCTION_H
#define NEWMESHREG_DISCRETEGROUPCOSTFUNCTION_H

#include "DiscreteCostFunction.h"

namespace newmeshreg {

class DiscreteGroupCostFunction : public NonLinearSRegDiscreteCostFunction {

public:
    DiscreteGroupCostFunction() = default;

    void set_parameters(myparam& p) override;
/*
    void set_relations(const boost::shared_ptr<RELATIONS>& CONTROL, const boost::shared_ptr<RELATIONS>& TARG){
        _controlrel = CONTROL;
        _targetrel = TARG;
        _sourcerel = _controlrel->invert_relationsTR( _CONTROLMESHES[0],_DATAMESHES[0]);
    }
*/
    void get_spacings();
    void set_meshes(const newresampler::Mesh& target, const newresampler::Mesh& source, const newresampler::Mesh& GRID, int num = 1) {
        _TEMPLATE = target;
        _ORIG = source;
        num_subjects = num;
        VERTICES_PER_SUBJ = GRID.nvertices();
        TRIPLETS_PER_SUBJ = GRID.ntriangles();
        MVD_LR = GRID.calculate_MeanVD();
        for(int i = 0; i < num; ++i) {
            _DATAMESHES.push_back(source);
            _CONTROLMESHES.push_back(GRID);
        }
    }

    //---INITIALISATION---//
    void initialize(int numNodes, int numLabels, int numPairs, int numTriplets); // quartets not used yet so no code for them below
    void define_template_patches();
    void resample_to_template();

    //---Updates---//
    void reset_source(const newresampler::Mesh& source, int num = 0) override { _DATAMESHES[num] = source; }
    virtual void reset_CPgrid(const newresampler::Mesh& grid, int num = 0) { _CONTROLMESHES[num] = grid; }

    double computeTripletCost(int triplet, int labelA, int labelB, int labelC) override;
    double computePairwiseCost(int pair, int labelA, int labelB) override;

    void resample_patches();
    void resampler_worker_function(int, int, const std::vector<bool>&);
    std::map<int,float> resample_onto_template(int, int, const newresampler::Point&, const std::vector<int>&);

private:
    std::vector<newresampler::Mesh> _DATAMESHES; // TARGET MESH
    std::vector<newresampler::Mesh> _CONTROLMESHES; // TARGET MESH
    newresampler::Mesh _TEMPLATE;

    std::vector<std::vector<int>> TEMPLATEPTS;
    std::vector<NEWMAT::ColumnVector> SPACINGS;

    //---DATA---//
    std::vector<NEWMAT::Matrix> RESAMPLEDDATA;
    std::vector<std::map<int,float>> PATCHDATA;

    double MVD_LR = 0.0;

    float _sigma = 0.0;
    float _lambdapairs = 1.0;
    int num_subjects = 0;
    int TRIPLETS_PER_SUBJ = 0;
    int VERTICES_PER_SUBJ = 0;

    bool _setpairs = false;
};

} //namespace newmeshreg

#endif //NEWMESHREG_DISCRETEGROUPCOSTFUNCTION_H
