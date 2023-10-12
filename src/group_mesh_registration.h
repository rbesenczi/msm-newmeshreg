#ifndef NEWMESHREG_GROUPMESHREG_H
#define NEWMESHREG_GROUPMESHREG_H

#include "mesh_registration.h"
#include "DiscreteGroupModel.h"
#include "newresampler/mesh.h"

namespace newmeshreg {

class Group_Mesh_registration : public Mesh_registration {

    std::vector<newresampler::Mesh> ALL_SPH_REG;
    newresampler::Mesh templ;
    int num_subjects = 1;

public:
    void initialize_level(int current_lvl) override;
    void evaluate() override;
    void transform(const std::string& filename) override;
    void run_discrete_opt(std::vector<newresampler::Mesh>&);
    void save_transformed_data(const std::string& filename) override;

    inline void set_inputs(const std::string& s) {
        std::vector<std::string> meshlist = read_ascii_list(s);
        num_subjects = meshlist.size();
        MESHES.clear();
        MESHES.resize(num_subjects);
        for (int subject = 0; subject < num_subjects; ++subject) {
            if(_verbose) std::cout << "Mesh #" << subject << " is " << meshlist[subject] << std::endl;
            MESHES[subject].load(meshlist[subject]);
            recentre(MESHES[subject]);
            true_rescale(MESHES[subject], RAD);
        }
    }

    inline void set_data_list(const std::string& s) {
        DATAlist = read_ascii_list(s);
        if (_verbose)
            for (int subject = 0; subject < DATAlist.size(); ++subject)
                std::cout << "Data #" << subject << " is " << DATAlist[subject] << std::endl;
    }

    inline void set_template(const std::string &M) {
        if(_verbose) std::cout << "Template is " << M << std::endl;
        templ.load(M);
        recentre(templ);
        true_rescale(templ,RAD);
    }

    inline void saveSPH_reg(const std::string& filename) const override {
        for(int subject = 0; subject < num_subjects; ++subject)
            ALL_SPH_REG[subject].save(filename + "sphere-" + std::to_string(subject) + ".LR.reg" + _surfformat);
    }
};

} //namespace newmeshreg

#endif //NEWMESHREG_GROUPMESHREG_H
