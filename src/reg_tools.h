#ifndef NEWMESHREG_REG_TOOLS_H
#define NEWMESHREG_REG_TOOLS_H

#include "newresampler/octree.h"
#include "miscmaths/histogram.h"

namespace newmeshreg {

class Neighbourhood {

    std::vector<std::vector<int>> neighbours;
    double angsep = 0.0;

public:
    Neighbourhood() = default;
    void update(const newresampler::Mesh& source, const newresampler::Mesh& target, double ang);

    inline unsigned long nrows(int i) const { return neighbours.at(i).size(); }
    inline int operator()(int i, int j) const { return neighbours.at(i).at(j); }
    std::vector<int>& at(int i) { return neighbours.at(i); }
    const std::vector<int>& at(int i) const { return neighbours.at(i); }
};

bool get_all_neighbours(int index, std::vector<int>& N, const newresampler::Point& point, int n,
                        const newresampler::Mesh& REF, std::shared_ptr<Neighbourhood>& _rel, MISCMATHS::SpMat<int>& found);

//---UNFOLD---//
void computeNormal2EdgeOfTriangle(const newresampler::Point& v0, const newresampler::Point& v1, const newresampler::Point& v2, newresampler::Point& norm2edge);
newresampler::Point computeGradientOfBarycentricTriangle(const newresampler::Point& v0, const newresampler::Point& v1, const newresampler::Point& v2);
newresampler::Point spatialgradient(int index, const newresampler::Mesh& SOURCE);
bool check_for_intersections(int ind, newresampler::Mesh &IN);
void unfold(newresampler::Mesh& SOURCE);

//---TANGS---//
newresampler::Tangs calculate_tangs(int ind, const newresampler::Mesh& SPH_in);
newresampler::Tangs calculate_tri(int ind, const newresampler::Mesh& SPH_in);
newresampler::Tangs calculate_tri(const newresampler::Point& a);

//---STRAINS---//
NEWMAT::Matrix get_coordinate_transformation(double dNdT1,double dNdT2, NEWMAT::ColumnVector& Norm);
NEWMAT::ColumnVector calculate_strains(int index, const std::vector<int>& kept, const newresampler::Mesh& orig, const newresampler::Mesh& final, const std::shared_ptr<NEWMAT::Matrix>& PrincipalStretches);
newresampler::Mesh calculate_strains(double fit_radius, const newresampler::Mesh& orig, const newresampler::Mesh& final, const std::shared_ptr<NEWMAT::Matrix>& PrincipalStretches = std::shared_ptr<NEWMAT::Matrix>());
double triangle_strain(const NEWMAT::Matrix& AA, const NEWMAT::Matrix & BB, double MU, double KAPPA, const std::shared_ptr<NEWMAT::ColumnVector>& strains, double k_exp);
double calculate_triangular_strain(int index, const newresampler::Mesh& ORIG, const newresampler::Mesh& FINAL, double mu, double kappa, const std::shared_ptr<NEWMAT::ColumnVector>& indexSTRAINS = std::shared_ptr<NEWMAT::ColumnVector>(), double k_exp = 2.0);
newresampler::Mesh calculate_triangular_strains(const newresampler::Mesh& ORIG, const newresampler::Mesh& FINAL, double MU, double KAPPA);
double calculate_triangular_strain(const newresampler::Triangle& ORIG_tr, const newresampler::Triangle& FINAL_tr, double mu, double kappa, const std::shared_ptr<NEWMAT::ColumnVector>& indexSTRAINS = std::shared_ptr<NEWMAT::ColumnVector>(), double k_exp = 2.0);

//---HIST NORMALISATION---//
void get_range(int dim, const MISCMATHS::BFMatrix& M, const NEWMAT::ColumnVector& excluded, double& min, double& max);
void set_range(int dim, MISCMATHS::BFMatrix& M, const NEWMAT::ColumnVector& excluded, double& min, double& max);
void multivariate_histogram_normalization(MISCMATHS::BFMatrix& IN, MISCMATHS::BFMatrix& REF, const std::shared_ptr<newresampler::Mesh>& EXCL_IN, const std::shared_ptr<newresampler::Mesh>& EXCL_REF, bool rescale = false);

//---READ DATA---//
void set_data(const std::string& dataname, std::shared_ptr<MISCMATHS::BFMatrix>& BF, newresampler::Mesh& M, bool issparse = false);

} // namespace newmeshreg

#endif //NEWMESHREG_REG_TOOLS_H
