#include "group_mesh_registration.h"

namespace newmeshreg {

void Group_Mesh_registration::initialize_level(int current_lvl) {

    check();
    if(cost[current_lvl] == "RIGID" || cost[current_lvl] == "AFFINE")
        throw MeshregException("AFFINE/RIGID registration is not supported in groupwise mode yet.");

    const std::vector<double> sigma(MESHES.size(), _sigma_in[current_lvl]);

    FEAT = std::make_shared<featurespace>(DATAlist);
    FEAT->set_smoothing_parameters(sigma);
    FEAT->set_cutthreshold(_threshold);
    FEAT->varnorm(_varnorm);
    FEAT->intensitynormalize(_IN, _cut);
    FEAT->is_sparse(_issparse);
    FEAT->set_nthreads(_numthreads);
    SPH_orig = FEAT->initialize(_genesis[current_lvl], MESHES, _exclude);
    if(FEAT->get_dim() > 1)
        throw MeshregException("Multivariate registration is not supported in groupwise mode yet.");
    PARAMETERS.insert(parameterPair("multivariate", false));

    newresampler::Mesh control = newresampler::make_mesh_from_icosa(_gridres[current_lvl]);
    newresampler::recentre(control);
    newresampler::true_rescale(control, RAD);

    model = std::make_shared<DiscreteGroupModel>(PARAMETERS);
    if(_debug) model->set_debug();
    model->set_featurespace(FEAT);
    model->set_meshspace(SPH_orig, SPH_orig, (int)MESHES.size());
    model->Initialize(control);
}

void Group_Mesh_registration::evaluate() {

    newresampler::Mesh oldreg;
    for(int i = 0; i < MESHES.size(); ++i)
    {
        if(level == 1)
            ALL_SPH_REG.push_back(project_CPgrid(SPH_orig, oldreg));
        else
        {
            oldreg = ALL_SPH_REG[i];
            ALL_SPH_REG[i] = project_CPgrid(SPH_orig, oldreg, i);
        }
    }

    run_discrete_opt(ALL_SPH_REG);

    if(_verbose) std::cout << "Exit main algorithm." << std::endl;
}

void Group_Mesh_registration::run_discrete_opt(std::vector<newresampler::Mesh>& meshes) {

    int iter = 1;
    double energy = 0.0, newenergy = 0.0;
    newresampler::Mesh transformed_controlgrid;
    std::vector<newresampler::Mesh> controlgrid;

    while(iter <= boost::get<int>(PARAMETERS.find("iters")->second))
    {
        controlgrid.clear();
        for(int i = 0; i < meshes.size(); ++i)
        {
            model->reset_meshspace(meshes[i], i);
            controlgrid.push_back(model->get_CPgrid(i));
        }
        model->setupCostFunction();
#ifdef HAS_HOCR
        newenergy = Fusion::optimize(model, _verbose, _numthreads);
#else
        throw MeshregException("Groupwise mode is only supported in the HOCR version of MSM.");
#endif
        if(iter > 3 && energy-newenergy < 0.001)
        {
            if (_verbose)
                std::cout << iter << " level has converged.\n"
                          << "newenergy " << newenergy << "\tenergy " << energy
                          << "\tenergy-newenergy " << energy - newenergy << std::endl;

            break;
        }

        model->applyLabeling();

        for(int i = 0; i < meshes.size(); ++i)
        {
            transformed_controlgrid = model->get_CPgrid(i);
            newresampler::barycentric_mesh_interpolation(meshes[i], controlgrid[i], transformed_controlgrid, _numthreads);
            unfold(transformed_controlgrid);
            model->reset_CPgrid(transformed_controlgrid, i);
            unfold(meshes[i]);
        }
        energy = newenergy;
        iter++;
    }
}

void Group_Mesh_registration::transform(const std::string &filename) {

    for(int i = 0; i < MESHES.size(); ++i)
    {
        newresampler::barycentric_mesh_interpolation(MESHES[i], SPH_orig, ALL_SPH_REG[i]);
        MESHES[i].save(filename + "sphere-" + std::to_string(i) + ".reg" + _surfformat);
    }
}

void Group_Mesh_registration::save_transformed_data(const std::string &filename) {
    for(int i = 0; i < MESHES.size(); ++i)
    {
        std::shared_ptr<MISCMATHS::BFMatrix> data;
        set_data(DATAlist[i], data, MESHES[i]);
        newresampler::metric_resample(MESHES[i], templ, _numthreads).save(filename + "transformed_and_reprojected-" + std::to_string(i) + _dataformat);
    }
}

} //namespace newmeshreg
