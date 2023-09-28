#ifndef NEWMESHREG_GROUPMESHREG_H
#define NEWMESHREG_GROUPMESHREG_H

#include "mesh_registration.h"
#include "DiscreteGroupModel.h"
#include "newresampler/mesh.h"

namespace newmeshreg {

class Group_Mesh_registration : public Mesh_registration {

    std::vector<newresampler::Mesh> ALL_SPH_REG;
    newresampler::Mesh templ;

public:
    void initialize_level(int current_lvl) override;
    void evaluate() override;
    void transform(const std::string& filename) override;
    void run_discrete_opt(std::vector<newresampler::Mesh>&);
    void save_transformed_data(const std::string& filename) override;

    inline void set_inputs(const std::string& s) {
        std::vector<std::string> meshlist = read_ascii_list(s);
        MESHES.clear();
        MESHES.resize(meshlist.size());
        for (int i = 0; i < meshlist.size(); ++i)
        {
            if(_verbose) std::cout << "Mesh #" << i << " is " << meshlist[i] << std::endl;
            newresampler::Mesh tmp;
            tmp.load(meshlist[i]);
            MESHES[i] = tmp;
        }
    }

    inline void set_data_list(const std::string& s) {
        DATAlist = read_ascii_list(s);
        if (_verbose)
            for (int i = 0; i < DATAlist.size(); ++i)
                std::cout << "Data #" << i << " is " << DATAlist[i] << std::endl;
    }

    inline void set_template(const std::string &M) {
        if(_verbose) std::cout << "Template is " << M << std::endl;
        templ.load(M);
        recentre(templ);
        true_rescale(templ,RAD);
    }

    inline void saveSPH_reg(const std::string& filename) const override {
        for(int i = 0; i < ALL_SPH_REG.size(); ++i)
            ALL_SPH_REG[i].save(filename + "sphere-" + std::to_string(i) + ".LR.reg" + _surfformat);
    }
};

} //namespace newmeshreg

#endif //NEWMESHREG_GROUPMESHREG_H
