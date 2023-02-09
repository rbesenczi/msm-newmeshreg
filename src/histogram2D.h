#ifndef NEWMESHREG_HISTOGRAM2D_H
#define NEWMESHREG_HISTOGRAM2D_H

#include "armawrap/newmat.h"
#include "miscmaths/miscmaths.h"

namespace newmeshreg {

class Histogram2D {

    int m_nbinsx;
    int m_nbinsy;
    int _nsamps;
    double m_min_x;
    double m_min_y;
    double m_max_x;
    double m_max_y;
    double m_width_x;
    double m_width_y;
    NEWMAT::Matrix _bins;
    NEWMAT::Matrix _weights;

public:
    inline double operator()(int i, int j){ return _bins(i,j); }

    void Initialize(int, int, double, double, double, double);

    // this can be used to remove bin matrix if memory overheads are a consideration
    // (i.e. when using a large number of separate histograms as in Costfunction.cc)
    void Reset() { _bins.ReSize(0,0); _nsamps = 0; }
    void Zero();
    void ReSize(int m, int n) { _bins.ReSize(m,n); _bins = 0; _nsamps = 0; }

    // add count for sample x,y
    void AddSample(double, double);
    void AddSample(double, double, double); // for weighted NMI

    // delete count for sample x,y
    void DeleteSample(double, double);

    // Entropy calculations
    double MarginalEntropyX();
    double MarginalEntropyY();
    double JointEntropy();

    // Normalised mutual Information
    double normalisedmutualinformation();
};

} //namespace newmeshreg

#endif //NEWMESHREG_HISTOGRAM2D_H
