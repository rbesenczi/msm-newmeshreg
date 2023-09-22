#ifndef NEWMESHREG_REG_TOOLS_H
#define NEWMESHREG_REG_TOOLS_H

#include "newresampler/octree.h"
#include "miscmaths/histogram.h"
#include <omp.h>
#include <boost/variant/variant.hpp>
#include <boost/variant/get.hpp>

namespace newmeshreg {

typedef std::map<std::string, boost::variant<int, std::string, double, float, bool>> myparam;
typedef std::pair<std::string, boost::variant<int, std::string, double, float, bool>> parameterPair;

class MeshregException : public std::exception {

public:
    const char* errmesg;

    explicit MeshregException(const char* msg) : errmesg(msg) {}

    const char* what() const noexcept override;
};

class Neighbourhood {

    std::vector<std::vector<int>> neighbours;
    double angsep = 0.0;

public:
    Neighbourhood() = default;
    void update(const newresampler::Mesh& source, const newresampler::Mesh& target, double ang, int numthreads = 1);

    inline unsigned long nrows(int i) const { return neighbours[i].size(); }
    inline int operator()(int i, int j) const { return neighbours[i][j]; }
    std::vector<int>& at(int i) { return neighbours[i]; }
    const std::vector<int>& at(int i) const { return neighbours[i]; }
};

bool get_all_neighbours(int index, std::vector<int>& N, const newresampler::Point& point, int n,
                        const newresampler::Mesh& REF, std::shared_ptr<Neighbourhood>& _rel, MISCMATHS::SpMat<int>& found);

//---UNFOLD---//
void computeNormal2EdgeOfTriangle(const newresampler::Point& v0, const newresampler::Point& v1, const newresampler::Point& v2, newresampler::Point& norm2edge);
newresampler::Point computeGradientOfBarycentricTriangle(const newresampler::Point& v0, const newresampler::Point& v1, const newresampler::Point& v2);
newresampler::Point spatialgradient(int index, const newresampler::Mesh& SOURCE);
bool check_for_intersections(int ind, newresampler::Mesh &IN);
void unfold(newresampler::Mesh& SOURCE, bool verbosity = false);

//---TANGS---//
newresampler::Tangs calculate_tangs(int ind, const newresampler::Mesh& SPH_in);
newresampler::Tangs calculate_tri(int ind, const newresampler::Mesh& SPH_in);
newresampler::Tangs calculate_tri(const newresampler::Point& a);

//---STRAINS---//
NEWMAT::Matrix get_coordinate_transformation(double dNdT1,double dNdT2, NEWMAT::ColumnVector& Norm);
NEWMAT::ColumnVector calculate_strains(int index, const std::vector<int>& kept, const newresampler::Mesh& orig, const newresampler::Mesh& final, const std::shared_ptr<NEWMAT::Matrix>& PrincipalStretches);
newresampler::Mesh calculate_strains(double fit_radius, const newresampler::Mesh& orig, const newresampler::Mesh& final, int numthreads = 1, const std::shared_ptr<NEWMAT::Matrix>& PrincipalStretches = std::shared_ptr<NEWMAT::Matrix>());
double triangle_strain(const NEWMAT::Matrix& AA, const NEWMAT::Matrix & BB, double MU, double KAPPA, const std::shared_ptr<NEWMAT::ColumnVector>& strains, double k_exp);
double calculate_triangular_strain(int index, const newresampler::Mesh& ORIG, const newresampler::Mesh& FINAL, double mu, double kappa, const std::shared_ptr<NEWMAT::ColumnVector>& indexSTRAINS = std::shared_ptr<NEWMAT::ColumnVector>(), double k_exp = 2.0);
double calculate_triangular_strain(const newresampler::Triangle& ORIG_tr, const newresampler::Triangle& FINAL_tr, double mu, double kappa, const std::shared_ptr<NEWMAT::ColumnVector>& indexSTRAINS = std::shared_ptr<NEWMAT::ColumnVector>(), double k_exp = 2.0);

//---HIST NORMALISATION---//
void multivariate_histogram_normalization(MISCMATHS::BFMatrix& IN, MISCMATHS::BFMatrix& REF, const std::shared_ptr<newresampler::Mesh>& EXCL_IN, const std::shared_ptr<newresampler::Mesh>& EXCL_REF, int nthreads = 1);

//---READ DATA---//
void set_data(const std::string& dataname, std::shared_ptr<MISCMATHS::BFMatrix>& BF, newresampler::Mesh& M, bool issparse = false);

template<typename Iterator> inline bool next_combination(const Iterator first, Iterator k, const Iterator last);

} // namespace newmeshreg

#endif //NEWMESHREG_REG_TOOLS_H
