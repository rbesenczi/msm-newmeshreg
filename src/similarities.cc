#include "similarities.h"

using namespace std;

namespace newmeshreg {

void sparsesimkernel::Initialize(int simval) {

    _sim = simval;

    if (m_A == nullptr) throw newmeshreg::MeshregException("SIMILARITIES:: Connectivity matrices have not been initliased");

    if(_sim == 1 || _sim == 2)
    {
        _rmeanA = meanvector(*m_A);
        if(m_A->Nrows() == 1) _meanA = _rmeanA(1);
        if(m_B != nullptr)
            _rmeanB = meanvector(*m_B);
        else
            _rmeanB = _rmeanA;

        if(m_B->Nrows() == 1) _meanB = _rmeanB(1);
    }

    // for NMI initialise histogram
    if(_sim == 3)
    {
        calc_range(m_A,maxx,minx);
        if(m_B == nullptr)
        {
            miny = minx;
            maxy = maxx;
        }
        else
            calc_range(m_B,maxy,miny);

        if (m_A->Nrows() < 256)
            hist.Initialize(m_A->Nrows() - 1, m_A->Nrows() - 1, maxx, maxy, minx, miny);
        else
            hist.Initialize(256, 256, maxx, maxy, minx, miny);
    }
}

double sparsesimkernel::get_sim_for_min(const vector<double>& input, const vector<double>& reference, const vector<double>& weights){

    double val = 0.0;
    switch (_sim)
    {
        case 1:
            val = SSD(input,reference, weights);
            break;
        case 2:
            val = corr(input,reference, weights);
            val = 1 - ( 1 + val) / 2; // scale between 0 and 1
            break;
        case 3:
            val = NMI(input,reference);
            break;
    }
    return val;
}

void sparsesimkernel::calculate_sim_column_nbh(int ind) {

    for (int j = 0; j < nbh->nrows(ind); j++)
        if ((*nbh)(ind, j) != 0)
            switch (_sim)
            {
                case 1:
                    mp->Set((int) (*nbh)(ind, j) + 1, ind + 1, -SSD(ind + 1, (int) (*nbh)(ind, j) + 1));
                    break;
                case 2:
                    mp->Set((int) (*nbh)(ind, j) + 1, ind + 1, corr(ind + 1, (int) (*nbh)(ind, j) + 1));
                    break;
                case 3:
                    mp->Set((int) (*nbh)(ind, j) + 1, ind + 1, NMI(ind + 1, (int) (*nbh)(ind, j) + 1));
                    break;
            }
}

// for the case where we are working on vector data and we have no knowledge of the full data m_A (currently used in discrete opt)
void sparsesimkernel::Initialize(int simval, const vector<double>& inputdata, const vector<double>& refdata, const vector<double>& weights) {

    _sim = simval;

    if(_sim == 3)
    {
        calc_range(inputdata,maxx,minx);
        calc_range(refdata,maxy,miny);
        if (inputdata.size() < 256)
            hist.Initialize(inputdata.size() - 1, inputdata.size() - 1, maxx, maxy, minx, miny);
        else
            hist.Initialize(256, 256, maxx, maxy, minx, miny);
    }
    else
    {
        _meanA = 0;
        _meanB = 0;
        double sum = 0;

        for(unsigned int i = 0; i < inputdata.size(); i++)
        {
            double varwght = 0;
            if(weights.size() == inputdata.size())
                varwght = weights[i];
            else
                varwght = 1;

            _meanA += varwght * inputdata[i];
            _meanB += varwght * refdata[i];
            sum += varwght;
        }
        if(sum > 0)
        {
            _meanA/=sum;
            _meanB/=sum;
        }
    }
}

NEWMAT::RowVector sparsesimkernel::meanvector(const MISCMATHS::BFMatrix& fdt_matrix) {

    NEWMAT::RowVector mean(fdt_matrix.Ncols());
    mean = 0;

    if (fdt_matrix.Nrows() == 1)
    {
        double sum = 0;
        for (unsigned int i = 0; i < fdt_matrix.Ncols(); i++)
            sum = sum + fdt_matrix.Peek(1, i + 1);
        for (unsigned int i = 0; i < fdt_matrix.Ncols(); i++)
            mean(i + 1) = sum / fdt_matrix.Ncols();
    }
    else
        for (unsigned int i = 0; i < fdt_matrix.Ncols(); i++)
        {
            double sum = 0;
            for (unsigned int j = 0; j < fdt_matrix.Nrows(); j++)
                sum = sum + fdt_matrix.Peek(j + 1, i + 1);
            mean(i + 1) = sum / fdt_matrix.Nrows();
        }

    return mean;
}

void sparsesimkernel::calc_range(const std::shared_ptr<MISCMATHS::BFMatrix>& mat, double& max, double& min) {

    for (unsigned int i = 1; i <= mat->Nrows(); i++)
        for (unsigned int j = 1; j <= mat->Ncols(); j++)
        {
            if (mat->Peek(i, j) > max) max = mat->Peek(i, j);
            if (mat->Peek(i, j) < min) min = mat->Peek(i, j);
        }
}

void sparsesimkernel::calc_range(const vector<double> &mat, double &max, double &min) {
    // probably change to STL's min/max search
    max = 0;
    min = 1e7;
    for (double i : mat)
    {
        if(i > max) max = i;
        if(i < min) min = i;
    }
}

double sparsesimkernel::corr(int i,int j){

    int num = 0, numB = 0;
    double prod = 0.0, varA = 0.0, varB = 0.0;

    if(m_B == nullptr) m_B = m_A;

    if(_rmeanA.Ncols() == 0)
    {
        _rmeanA=meanvector(*m_A);
        if(_rmeanB.Ncols()==0)
        {
            if(m_B == nullptr)
                _rmeanB=_rmeanA;
            else
                _rmeanB=meanvector(*m_B);
        }
    }

    std::shared_ptr<MISCMATHS::SparseBFMatrix<double>> ptr = std::dynamic_pointer_cast<MISCMATHS::SparseBFMatrix<double> >(m_A);

    double Bzerooffset = (0.0-_rmeanB(j));  // result for all zero values of sparse mat
    double Azerooffset = (0.0-_rmeanA(i));

    for (MISCMATHS::BFMatrixColumnIterator it = m_A->begin(i); it != m_A->end(i); it++)
    {
        prod += (*it-_rmeanA(i))*(m_B->Peek(it.Row(),j)-_rmeanB(j));
        varA += (*it-_rmeanA(i))*(*it-_rmeanA(i));
        varB += (m_B->Peek(it.Row(),j)-_rmeanB(j))*(m_B->Peek(it.Row(),j)-_rmeanB(j));
        num++;
    }

    if(ptr)
    { //if sparse run all for all rows where A had zero values

        varA += (m_A->Nrows()-num)*Azerooffset*Azerooffset; // for rows where A has no values
        for (MISCMATHS::BFMatrixColumnIterator it = m_B->begin(j); it != m_B->end(j); it++)
        {
            if(m_A->Peek(it.Row(),i) == 0)
            {
                prod +=Azerooffset*(*it-_rmeanB(j));
                num++;
                varB+=(*it-_rmeanB(j))*(*it-_rmeanB(j));
            }
            numB++;
        }
    varB += (m_A->Nrows()-numB)*Bzerooffset*Bzerooffset;    // for rows where B has no values
    prod += (2*m_A->Nrows()-num)*Azerooffset*Bzerooffset;  // for rows where A&B have no values
    }

    if (varA == 0.0 || varB == 0.0) return 0.0;
    else return prod/(sqrt(varA)*sqrt(varB));
}

double sparsesimkernel::SSD(int i, int j) {

    double prod = 0;
    if(m_B == nullptr)  m_B = m_A;

    std::shared_ptr<MISCMATHS::SparseBFMatrix<double>> ptr = std::dynamic_pointer_cast<MISCMATHS::SparseBFMatrix<double>>(m_A);

    for (MISCMATHS::BFMatrixColumnIterator it = m_A->begin(i); it != m_A->end(i); it++)
      prod += (*it-m_B->Peek(it.Row(),j))*(*it-m_B->Peek(it.Row(),j));

    if (ptr)
        for (MISCMATHS::BFMatrixColumnIterator it = m_B->begin(j); it != m_B->end(j); it++)
            if (m_A->Peek(it.Row(), i) == 0)
                prod += (*it) * (*it);

    return sqrt(prod)/m_A->Nrows();
}

/* NMI derived similarity kernel*/
double sparsesimkernel::NMI(int i, int j) {

    if(m_A->Nrows() < 256)
        hist.ReSize(m_A->Nrows()-1,m_A->Nrows()-1);
    else
        hist.ReSize(256,256);

    if(m_B == nullptr) m_B = m_A;

    for (unsigned int s = 1; s <= m_A->Nrows(); s++)
        hist.AddSample(m_A->Peek(s, i), m_B->Peek(s, j));

    return hist.normalisedmutualinformation();
}

//for vectors - weight function is optional
double sparsesimkernel::corr(const vector<double>& A,const vector<double>& B,const vector<double>& weights) {

    double prod = 0.0, varA = 0.0, varB = 0.0;
    double sum = 0.0;

    if(A.size() != B.size())
    {
        cout << A.size() << " " << B.size() << " " << weights.size() << endl;
        throw newmeshreg::MeshregException("SIMILARITIES:: correlation, data dimensions do not match");
    }

    Initialize(2,A,B, weights);

    for (unsigned int s = 0; s < A.size(); s++)
    {
        double varwght = 0;

        if (weights.size() == A.size()) varwght = weights[s];
        else varwght = 1;

        prod = prod + varwght * (A[s] - _meanA) * (B[s] - _meanB);
        varA = varA + varwght * (A[s] - _meanA) * (A[s] - _meanA);
        varB = varB + varwght * (B[s] - _meanB) * (B[s] - _meanB);

        sum += varwght;
    }

    if(sum > 0)
    {
      prod /= sum;
      varA /= sum;
      varB /= sum;
    }

    if (varA == 0.0 || varB == 0.0) return 0.0;
    else return prod / (sqrt(varA) * sqrt(varB));
}

double sparsesimkernel::SSD(const vector<double> &A, const vector<double> &B, const vector<double> &weights){

    double prod = 0.0, sum = 0.0;

    if (A.size() != B.size()) throw newmeshreg::MeshregException("SIMILARITIES:: SSD data dimensions do not match");

    Initialize(3,A,B,weights);

    for (unsigned int s = 0; s < A.size(); s++)
    {
        double varwght = 0;
        if(weights.size() == A.size()) varwght = weights[s];
        else varwght = 1;
        prod = prod + varwght * (A[s] - B[s])*(A[s] - B[s]);
        sum += varwght;
    }
    if(sum > 0.0) prod = prod/sum;
    return prod;
}

double sparsesimkernel::NMI(const vector<double> &A,const vector<double> &B,const vector<double> &weights){

    if (A.size() != B.size()) throw newmeshreg::MeshregException("SIMILARITIES:: NMI data dimensions do not match");

    Initialize(4,A,B,weights);

    if(A.size() < 256)
        hist.ReSize(A.size()-1,A.size()-1);
    else
        hist.ReSize(256,256);

    for (unsigned int s = 0; s < A.size(); s++)
    {
        double varwght = 0.0;
        if(weights.size() == A.size()) varwght = 1;
        else varwght = 1;
        hist.AddSample(varwght,A[s],B[s]);
    }

    return hist.normalisedmutualinformation();
}

} //namespace newmeshreg
