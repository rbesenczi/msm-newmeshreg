#ifndef NEWMESHREG_DISCRETEGROUPMODEL_H
#define NEWMESHREG_DISCRETEGROUPMODEL_H

#include "DiscreteModel.h"
#include "DiscreteGroupCostFunction.h"

namespace newmeshreg {

class DiscreteGroupModel : public NonLinearSRegDiscreteModel {

    std::vector<newresampler::Mesh> m_datameshes;
    std::vector<newresampler::Mesh> m_controlmeshes;
    newresampler::Mesh m_template;
    std::vector<std::map<int,double>> patch_data;
    std::vector<NEWMAT::ColumnVector> spacings;
    std::vector<NEWMAT::Matrix> rotated_meshes;

    int m_num_subjects = 0;
    int control_grid_size = 0;

public:
    DiscreteGroupModel() = default;
    explicit DiscreteGroupModel(myparam& p) {
        set_parameters(p);
        costfct = std::make_shared<DiscreteGroupCostFunction>();
        costfct->set_parameters(p);
    }

    void set_meshspace(const newresampler::Mesh& target, const newresampler::Mesh& source, int num = 1) override {
        m_template = target;
        m_datameshes.clear();
        m_datameshes.resize(num, source);
        m_num_subjects = num;
    }

    void reset_meshspace(const newresampler::Mesh& source, int num = 0) override {
        m_datameshes[num] = source;
        //costfct->reset_source(source, num);
    }

    void reset_CPgrid(const newresampler::Mesh& grid, int num = 0) override { m_controlmeshes[num] = grid; }

    void warp_CPgrid(newresampler::Mesh& start, newresampler::Mesh& end, int num = 0) override {
        newresampler::barycentric_mesh_interpolation(m_controlmeshes[num], start, end, _nthreads);
        unfold(m_controlmeshes[num], m_verbosity);
    }

    void applyLabeling() override { applyLabeling(labeling); }
    void applyLabeling(int* dlabels) override {
        #pragma omp parallel for num_threads(_nthreads)
        for (int subject = 0; subject < m_num_subjects; subject++)
            for (int vertex = 0; vertex < control_grid_size; vertex++)
                m_controlmeshes[subject].set_coord(vertex, m_ROT[subject * control_grid_size + vertex] *
                                                           m_labels[dlabels[subject * control_grid_size + vertex]]);
    }

    newresampler::Mesh get_CPgrid(int num = 0) override { return m_controlmeshes[num]; }

    void initialize_pairs();
    void estimate_pairs() override;
    void estimate_triplets() override;
    void get_patch_data();
    void get_rotated_meshes();

    void Initialize(const newresampler::Mesh& controlgrid) override;
    void get_rotations(std::vector<NEWMAT::Matrix>&) override;
    void setupCostFunction() override;
};

} //namespace newmeshreg

#endif //NEWMESHREG_DISCRETEGROUPMODEL_H
