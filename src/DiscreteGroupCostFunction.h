#ifndef NEWMESHREG_DISCRETEGROUPCOSTFUNCTION_H
#define NEWMESHREG_DISCRETEGROUPCOSTFUNCTION_H

#include "DiscreteCostFunction.h"

namespace newmeshreg {

class DiscreteGroupCostFunction : public NonLinearSRegDiscreteCostFunction {

public:
    void set_meshes(const newresampler::Mesh& target, const newresampler::Mesh& source, const newresampler::Mesh& GRID, int num = 1) {
        _TEMPLATE = target;
        _ORIG = source;
        num_subjects = num;
        VERTICES_PER_SUBJ = GRID.nvertices();
        TRIPLETS_PER_SUBJ = GRID.ntriangles();
        _DATAMESHES.resize(num, source);
        _CONTROLMESHES.resize(num,GRID);
    }
    //---INITIALISATION---//
    void initialize(int numNodes, int numLabels, int numPairs, int numTriplets) override;
    void set_group_spacings(const std::vector<NEWMAT::ColumnVector>& sp) override { SPACINGS = sp; }
    void set_targettree(const std::shared_ptr<newresampler::Octree>& tree) override { targettree = tree; }

    //---Updates---//
    void reset_source(const newresampler::Mesh& source, int num = 0) override { _DATAMESHES[num] = source; }
    virtual void reset_CPgrid(const newresampler::Mesh& grid, int num = 0) { _CONTROLMESHES[num] = grid; }

    //---Compute costs---//
    double computeTripletCost(int triplet, int labelA, int labelB, int labelC) override;
    double computePairwiseCost(int pair, int labelA, int labelB) override;
    void get_patch_data();
    void set_rotated_meshes(const std::vector<newresampler::Mesh>& MESHES, const std::vector<NEWMAT::Matrix>& data) override {
        _ROTATED_DATAMESHES = MESHES;
        resampled_data = data;
    }
    void set_patch_data(const std::vector<std::map<int,double>>& patches) override {
        patch_data = patches;
    }

private:
    std::vector<newresampler::Mesh> _DATAMESHES;
    std::vector<newresampler::Mesh> _ROTATED_DATAMESHES;
    std::vector<newresampler::Mesh> _CONTROLMESHES;
    newresampler::Mesh _TEMPLATE;
    std::vector<std::shared_ptr<newresampler::Octree>> datameshtrees;
    std::vector<NEWMAT::Matrix> resampled_data;

    std::vector<NEWMAT::ColumnVector> SPACINGS;

    std::vector<std::map<int,double>> patch_data;

    int num_subjects = 0;
    int TRIPLETS_PER_SUBJ = 0;
    int VERTICES_PER_SUBJ = 0;
};

} //namespace newmeshreg

#endif //NEWMESHREG_DISCRETEGROUPCOSTFUNCTION_H
