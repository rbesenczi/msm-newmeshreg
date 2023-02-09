#include "histogram2D.h"

namespace newmeshreg {

void Histogram2D::Initialize(int nbinsx, int nbinsy, double maxx, double maxy, double minx, double miny) {

    m_nbinsx  = nbinsx;
    m_nbinsy  = nbinsy;
    m_min_x   = minx;
    m_min_y   = miny;
    m_max_x   = maxx;
    m_max_y   = maxy;
    m_width_x = (m_max_x - m_min_x) / (double)m_nbinsx;
    m_width_y = (m_max_y - m_min_y) / (double)m_nbinsy;

    NEWMAT::Matrix newmat(m_nbinsx, m_nbinsy);

    _bins=newmat;
    _weights=newmat;
    _bins=0;
    _nsamps=0;
    _weights=0;
}

void Histogram2D::Zero() {

    for (int j = 1; j <= m_nbinsy; j++)
        for (int i = 1; i <= m_nbinsx; i++)
            _bins(i, j) = 0;

    _nsamps = 0;
}

void Histogram2D::AddSample(double x, double y) {

    if (x < m_min_x) return;
    if (x > m_max_x) return;
    if (y < m_min_y) return;
    if (y > m_max_y) return;

    int i = (int) MISCMATHS::round(m_nbinsx * (x - m_min_x - 0.5 * m_width_x) / (m_max_x - m_min_x));
    int j = (int) MISCMATHS::round(m_nbinsy * (y - m_min_y - 0.5 * m_width_y) / (m_max_y - m_min_y));

    if (i < 1) i = 1;
    if (j < 1) j = 1;
    if (i > m_nbinsx) i = m_nbinsx;
    if (j > m_nbinsy) j = m_nbinsy ;

    _bins(i,j) += 1;
    _weights(i,j) = 1;
    _nsamps += 1;
}

void Histogram2D::AddSample(double weight, double x, double y) {

    if (x < m_min_x) return;
    if (x > m_max_x) return;
    if (y < m_min_y) return;
    if (y > m_max_y) return;

    int i = (int) MISCMATHS::round(m_nbinsx * (x - m_min_x - 0.5*m_width_x) / (m_max_x - m_min_x));
    int j = (int) MISCMATHS::round(m_nbinsy * (y - m_min_y - 0.5*m_width_y) / (m_max_y - m_min_y));

    if (i < 1) i = 1;
    if (j < 1) j = 1;
    if (i > m_nbinsx) i = m_nbinsx;
    if (j > m_nbinsy) j = m_nbinsy ;

    _bins(i,j) += 1;
    _weights(i,j) = 1;// set as one for the time being until we can find correct formula weight;
    _nsamps += 1;
}

void Histogram2D::DeleteSample(double x, double y) {

    if (x < m_min_x) return;
    if (x > m_max_x) return;
    if (y < m_min_y) return;
    if (y > m_max_y) return;

    int i = (int) MISCMATHS::round(m_nbinsx * (x - m_min_x - 0.5*m_width_x) / (m_max_x - m_min_x));
    int j = (int) MISCMATHS::round(m_nbinsy * (y - m_min_y - 0.5*m_width_y) / (m_max_y - m_min_y));

    if (i < 0) i = 1;
    if (j < 0) j = 1;
    if (i > m_nbinsx) i = m_nbinsx;
    if (j > m_nbinsy) j = m_nbinsy;

    _bins(i,j) -= 1;
    _nsamps -= 1;
}

double Histogram2D::MarginalEntropyX() {

    NEWMAT::ColumnVector M(m_nbinsx);
    M=0;
    double E = 0.0;

    for (int i = 1; i <= m_nbinsx; i++)
        for (int j = 1; j <= m_nbinsy; j++)
            M(i) = M(i) + _bins(i, j) / (double) _nsamps;

    for (int i = 1; i <= m_nbinsx; i++)
        if (M(i) > 0)
            E -= M(i) * log(M(i));

    return E;
}

double Histogram2D::MarginalEntropyY() {

    NEWMAT::ColumnVector M(m_nbinsy);
    M=0;
    double E = 0.0;

    for (int j = 1; j <= m_nbinsy; j++)
        for (int i = 1; i <= m_nbinsx; i++)
            M(j) = M(j) + _bins(i, j) / (double) _nsamps;

    for (int i = 1; i <= m_nbinsx; i++)
        if (M(i) > 0)
            E -= M(i) * log(M(i));

    return E;
}

double Histogram2D::JointEntropy() {

    double E = 0.0;
    NEWMAT::Matrix JP(m_nbinsx,m_nbinsy);

    for (int i = 1; i <= m_nbinsy; i++)
        for (int j = 1; j <= m_nbinsx; j++)
            JP(i, j) = _bins(i, j) / (double) _nsamps;

    for (int i = 1; i <= m_nbinsy; i++)
        for (int j = 1; j <= m_nbinsx; j++)
            if (JP(i, j) > 0)
                E -= _weights(i, j) * JP(i, j) * log(JP(i, j));

    return E;
}

double Histogram2D::normalisedmutualinformation() {

    if(this->JointEntropy() > 1e-5)
        return (this->MarginalEntropyX() + this->MarginalEntropyY()) / this->JointEntropy();
    else return 0.0;
}

} //namespace newmeshreg
