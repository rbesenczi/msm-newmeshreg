#ifndef NEWMESHREG_DISCRETEGROUPMODEL_H
#define NEWMESHREG_DISCRETEGROUPMODEL_H

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

    int m_num_subjects = 0;
    int control_grid_size = 0;

public:
    DiscreteGroupModel() = default;
    explicit DiscreteGroupModel(myparam& p) {
        set_parameters(p);
        costfct = std::shared_ptr<NonLinearSRegDiscreteCostFunction>(new DiscreteGroupCostFunction());
        /*
         * relations
         * m_inputrel=boost::shared_ptr<RELATIONS >(new RELATIONS());
         * m_cp_neighbourhood=boost::shared_ptr<RELATIONS >(new RELATIONS ());
         */
        costfct->set_parameters(p);
    }

    void set_meshspace(const newresampler::Mesh& target, const newresampler::Mesh& source, int num = 1) override {
        m_template = target;
        m_datameshes.clear();
        m_num_subjects = num;
        for(int i = 0; i < num; ++i)
            m_datameshes.push_back(source);
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
    void applyLabeling(int* discreteLabeling) override;

    newresampler::Mesh get_CPgrid(int num = 0) override { return m_controlmeshes[num]; }

    void initialize_quartets();
    void initialize_pairs();

    void estimate_pairs() override;
    void estimate_triplets() override;
    void estimate_quartets();
    void estimate_combinations(int, int*);

    void Initialize(const newresampler::Mesh&) override;
    //void Initialize(){};
    void get_rotations(std::vector<NEWMAT::Matrix>&) override;
    void setupCostFunction() override;

    void get_between_subject_pairs();
    //void resample_to_template();
};

} //namespace newmeshreg

#endif //NEWMESHREG_DISCRETEGROUPMODEL_H
