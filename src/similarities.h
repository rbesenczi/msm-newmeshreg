#ifndef NEWMESHREG_SIMILARITIES_H
#define NEWMESHREG_SIMILARITIES_H

#include <miscmaths/bfmatrix.h>

#include <memory>
#include "meshregException.h"
#include "histogram2D.h"
#include "reg_tools.h"

namespace newmeshreg {

class sparsesimkernel{

public:
    sparsesimkernel(): maxx(0), maxy(0), _sim(1), _meanA(0.0), _meanB(0.0),
                       minx(std::numeric_limits<double>::max()),
                       miny(std::numeric_limits<double>::max()),
                       mp(std::make_shared<MISCMATHS::SpMat<double>>()){}

    //---INITIALISE---//
    void Initialize(int);

    //---ACCESS---//
    inline double Peek(unsigned int r, unsigned int c) const { return(mp->Peek(r,c)); }
    double get_sim_for_min(const std::vector<double>&, const std::vector<double>&, const std::vector<double>& weights = std::vector<double>());

    //---SET---//
    void set_simval(int val) { _sim = val; }
    void set_input(std::shared_ptr<MISCMATHS::BFMatrix> in ){ m_A = in; }
    void set_reference(std::shared_ptr<MISCMATHS::BFMatrix> ref){ m_B = ref; }
    void set_neighbourhood(std::shared_ptr<Neighbourhood>& n) { nbh = n; }
    void Resize(unsigned int m, unsigned int n) { mp = std::make_shared<MISCMATHS::SpMat<double>>(m,n); }
    void calculate_sim_column_nbh(int);

private:
    // DATA
    std::shared_ptr<MISCMATHS::BFMatrix> m_A, m_B;
    std::shared_ptr<MISCMATHS::SpMat<double>> mp;
    std::shared_ptr<Neighbourhood> nbh;

    NEWMAT::RowVector _rmeanA, _rmeanB;

    double _meanA, _meanB, maxx, maxy, minx, miny;
    int _sim;
    Histogram2D hist; // for NMI

    void Initialize(int, const std::vector<double>&, const std::vector<double>&, const std::vector<double>& weights = std::vector<double>()); // weighted version
    //---CALC---//
    NEWMAT::RowVector meanvector(const MISCMATHS::BFMatrix &); // for correlation measure
    void calc_range(const std::shared_ptr<MISCMATHS::BFMatrix>&, double&, double&); // for histogram based measures
    void calc_range(const std::vector<double>&, double&, double&);
    double corr(int, int);
    double SSD(int, int);
    double NMI(int, int);
    double corr(const std::vector<double>&, const std::vector<double>&, const std::vector<double>& weights = std::vector<double>());
    double SSD(const std::vector<double>&,const std::vector<double>&,const std::vector<double>& weights = std::vector<double>());
    double NMI(const std::vector<double>&, const std::vector<double>&, const std::vector<double>& weights = std::vector<double>());
};

} //namespace newmeshreg

#endif //NEWMESHREG_SIMILARITIES_H
