#ifndef NEWMESHREG_DISCRETEGROUPMODEL_H
#define NEWMESHREG_DISCRETEGROUPMODEL_H

#include <memory>

#include "DiscreteModel.h"
#include "DiscreteGroupCostFunction.h"

namespace newmeshreg {

class DiscreteGroupModel : public NonLinearSRegDiscreteModel {

    std::vector<newresampler::Mesh> m_datameshes;
    std::vector<newresampler::Mesh> m_controlmeshes;
    newresampler::Mesh m_template;
    newresampler::Mesh m_template_LR;

    std::vector<std::vector<std::vector<int>>> between_subject_pairs;
    std::vector<int> subjects;
    std::vector<std::shared_ptr<newresampler::Octree>> datameshtrees;

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
        costfct->reset_source(source, num);
    }

    void reset_CPgrid(const newresampler::Mesh& grid, int num = 0) override { m_controlmeshes[num] = grid; }

    void warp_CPgrid(newresampler::Mesh& start, newresampler::Mesh& end, int num = 0) override {
        newresampler::barycentric_mesh_interpolation(m_controlmeshes[num], start, end);
        unfold(m_controlmeshes[num]);
    }

    void applyLabeling() override { applyLabeling(labeling); }
    void applyLabeling(int* dlabels) override {
        for (int n = 0; n < m_num_subjects; n++)
            for (int i = 0; i < m_controlmeshes[n].nvertices(); i++)
                m_controlmeshes[n].set_coord(i, m_ROT[i + n * m_controlmeshes[n].nvertices()] *
                                                m_labels[dlabels[i + n * m_controlmeshes[n].nvertices()]]);
    }

    newresampler::Mesh get_CPgrid(int num = 0) override { return m_controlmeshes[num]; }

    void initialize_pairs();
    void estimate_pairs() override;
    void estimate_triplets() override;

    void Initialize(const newresampler::Mesh& controlgrid) override;
    void get_rotations(std::vector<NEWMAT::Matrix>&) override;
    void setupCostFunction() override;

    void get_between_subject_pairs();
};

} //namespace newmeshreg

#endif //NEWMESHREG_DISCRETEGROUPMODEL_H
