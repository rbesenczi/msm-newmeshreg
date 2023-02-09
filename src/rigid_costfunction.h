#ifndef NEWMESHREG_RIGID_COSTFUNCTION_H
#define NEWMESHREG_RIGID_COSTFUNCTION_H

#include "armawrap/newmat.h"
#include "miscmaths/SpMat.h"
#include "featurespace.h"
#include "similarities.h"
#include "newresampler/octree.h"
#include "reg_tools.h"

#ifdef HAS_HOCR
#include "Fusion.h"
#else
#include <FastPD/FastPD.h>
#endif

namespace newmeshreg {

//struct Tangs { newresampler::Point e1, e2; };

class Rigid_cost_function {
// basic mesh cost function - uses the gradient of the similarity function, estimated using weighted least squares

    newresampler::Mesh TARGET; // TARGET MESH
    newresampler::Mesh SOURCE; // SOURCE MESH

    std::shared_ptr<Neighbourhood> nbh;
    std::shared_ptr<newresampler::Octree> targettree;
    NEWMAT::Matrix inweight; // exclusion/weighting mask for cost function masking
    NEWMAT::Matrix refweight;

    std::shared_ptr<featurespace> FEAT; // holds data
    sparsesimkernel sim; // similarity matrix
    std::vector<newresampler::Tangs> BASIS; // tangent vector bases for all vertices
    NEWMAT::ColumnVector REC;

    double MVD; // mean vertex distance
    NEWMAT::ColumnVector current_sim;
    double min_sigma;

    // user defined parameters
    int simmeasure;
    int iters; // total iterations
    float stepsize;
    float spacing;

    //std::string m_outdir;
    bool verbosity;

    //---SIMILARITY GRADIENT ESTIMATION---//
    // similarity gradient via weighted regression
    NEWMAT::ColumnVector WLS_simgradient(const newresampler::Tangs& tangs, int index, const std::vector<int>& querypoints);
    // prepares data for sim gradient calculation
    NEWMAT::ColumnVector Evaluate_SIMGradient(int i, const newresampler::Tangs& tangs);

    //---TRANSFORM AND EVALUATE---//
    void rotate_in_mesh(double a1, double a2, double a3);
    double rigid_cost_mesh(double dw1, double dw2, double dw3);

    //---MAKE UPDATES---//
    void update_similarity();

public:
    Rigid_cost_function(const newresampler::Mesh& target, const newresampler::Mesh& source,
                        std::shared_ptr<featurespace>& features);
    void set_parameters(myparam& PAR);

    //---INITIALIZE AND UPDATE---//
    void Initialize();
    void set_simmeasure(int simval) { simmeasure = simval; }
    void update_source(const newresampler::Mesh& M) { SOURCE = M; }

    //---EXECUTE---//
    newresampler::Mesh run();
};

//Tangs calculate(int ind, const newresampler::Mesh& SPH_in); // should move this to point or Triangle?
//void project_point(const newresampler::Point& vb, const Tangs& T, double& e1coord, double& e2coord); // should move this to point or Triangle?

} //namespace newmeshreg

#endif //NEWMESHREG_RIGID_COSTFUNCTION_H
