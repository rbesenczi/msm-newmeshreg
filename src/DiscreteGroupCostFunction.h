#ifndef NEWMESHREG_DISCRETEGROUPCOSTFUNCTION_H
#define NEWMESHREG_DISCRETEGROUPCOSTFUNCTION_H

#include "DiscreteCostFunction.h"

namespace newmeshreg {

class DiscreteGroupCostFunction : public NonLinearSRegDiscreteCostFunction {

public:
    void set_parameters(myparam& p) override;
    void set_meshes(const newresampler::Mesh& target, const newresampler::Mesh& source, const newresampler::Mesh& GRID, int num = 1) {
        _TEMPLATE = target;
        _ORIG = source;
        num_subjects = num;
        VERTICES_PER_SUBJ = GRID.nvertices();
        TRIPLETS_PER_SUBJ = GRID.ntriangles();
        MVD_LR = GRID.calculate_MeanVD();
        _DATAMESHES.resize(num, source);
        _CONTROLMESHES.resize(num,GRID);
    }
    //---INITIALISATION---//
    void initialize(int numNodes, int numLabels, int numPairs, int numTriplets) override;
    void resample_to_template();
    void get_spacings();
    void set_trees(const std::vector<std::shared_ptr<newresampler::Octree>>& trees) override { datameshtrees = trees; }

    //---Updates---//
    void reset_source(const newresampler::Mesh& source, int num = 0) override { _DATAMESHES[num] = source; }
    virtual void reset_CPgrid(const newresampler::Mesh& grid, int num = 0) { _CONTROLMESHES[num] = grid; }

    //---Compute costs---//
    double computeTripletCost(int triplet, int labelA, int labelB, int labelC) override;
    double computePairwiseCost(int pair, int labelA, int labelB) override;
    std::vector<double> get_patch_data(int, const NEWMAT::Matrix& rot);
    void get_source_data() override;

private:
    std::vector<newresampler::Mesh> _DATAMESHES; // TARGET MESH
    std::vector<newresampler::Mesh> _CONTROLMESHES; // TARGET MESH
    newresampler::Mesh _TEMPLATE;
    std::vector<std::shared_ptr<newresampler::Octree>> datameshtrees;

    std::vector<NEWMAT::ColumnVector> SPACINGS;

    std::vector<NEWMAT::Matrix> RESAMPLEDDATA;

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
