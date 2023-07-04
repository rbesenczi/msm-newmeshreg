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
    void initialize_level(int current_lvl);
    void evaluate();
    void transform(const std::string& filename);
    void run_discrete_opt(std::vector<newresampler::Mesh>&);
    void save_transformed_data(const std::string& filename);

    inline void saveSPH_reg(const std::string& filename) const {
        for(int i = 0; i < ALL_SPH_REG.size(); ++i)
            ALL_SPH_REG[i].save(filename + "sphere-" + std::to_string(i) + ".LR.reg" + _surfformat);
    }

    inline void set_data_list(const std::string& s) { DATAlist = read_ascii_list(s); }
    inline void set_template(const std::string &M) {templ.load(M); recentre(templ);  true_rescale(templ,RAD); }

    void set_inputs(const std::string& s) {
        std::vector<std::string> meshlist = read_ascii_list(s);
        newresampler::Mesh tmp;
        MESHES.clear();
        for (int i = 0; i < meshlist.size(); i++) {
            if (_verbose) std::cout << i << " " << meshlist[i] << std::endl;
            tmp.load(meshlist[i]);
            MESHES.push_back(tmp);
        }
    }
};

} //namespace newmeshreg

#endif //NEWMESHREG_GROUPMESHREG_H
