#ifndef NEWMESHREG_FEATURESPACE_H
#define NEWMESHREG_FEATURESPACE_H

#include <fstream>
#include <cstdio>
#include <memory>
#include <vector>

#include "newresampler/resampler.h"
#include "meshregException.h"
#include "reg_tools.h"

#define RAD 100.0
#define EPSILON 1.0E-8

namespace newmeshreg {

class featurespace {

public:
    featurespace() = default;
    featurespace(const std::string& datain, const std::string& dataref);
    featurespace(const std::string& datain, const std::vector<std::string>& datareflist);
    explicit featurespace(const std::vector<std::string>& datalist);
    ~featurespace() = default;

    //---INITIALISE---//
    newresampler::Mesh initialize(int, std::vector <newresampler::Mesh> &, bool);

    //---SET---//
    void set_smoothing_parameters(const std::vector<double>& s);
    void set_cutthreshold(std::vector<float>& thr) { _fthreshold = thr; }
    void logtransform(bool log) { _logtransform = log; }
    void varnorm(bool norm) { _varnorm = norm; }
    void is_sparse(bool sp) { _issparse = sp; }
    void intensitynormalize(bool norm, bool _exclcut) { _intensitynorm = norm; _cut = _exclcut; }
    void resamplingmethod(const std::string& method) { _resamplingmethod = method; }

    //---GET---//
    int get_dim() const { return DATA[0]->Nrows(); }
    double get_input_val(int i, int j) const { return DATA[0]->Peek(i, j); }
    double get_ref_val(int i, int j) const { return DATA[1]->Peek(i,j); }
    std::shared_ptr<newresampler::Mesh> get_input_excl() const { return EXCL[0]; }
    std::shared_ptr<newresampler::Mesh> get_reference_excl() const { return EXCL[1]; }
    std::shared_ptr<MISCMATHS::BFMatrix> get_input_data() const { return DATA[0]; }
    std::shared_ptr<MISCMATHS::BFMatrix> get_reference_data() const { return DATA[1]; }

private:
    newresampler::Mesh source;
    std::vector<std::shared_ptr<MISCMATHS::BFMatrix>> DATA; // holds generic BFMATRIX data which can be sparse or full matrices
    std::vector<std::shared_ptr<newresampler::Mesh>> EXCL;  // exclusion masks for binary weighting of the data during resampling

    std::vector<std::string> CMfile_in;  // path to data
    std::vector<double> _sigma_in;  // smoothing parameters for input and reference
    std::vector<float> _fthreshold;

    std::string inorig;
    std::string reforig;
    std::string _resamplingmethod;

    bool _logtransform = false;  // will log transform and normalise
    bool _intensitynorm = false; // will histogram match
    bool _scale = false;  // will rescale each feature to crudley match the distribution of the first in a multivariate distribution
    bool _issparse = false;  // notes that data is sparse
    bool _varnorm = false;  // performs online variance normalisation - maybe replace with non online version called during logtransformandnormalise()
    bool _cut = false;

    //---PROCESSING---//
    void variance_normalise(std::shared_ptr<MISCMATHS::BFMatrix>&, std::shared_ptr<newresampler::Mesh>&);
    void log10_transform_and_normalise(MISCMATHS::BFMatrix& data);
};

} //namespace newmeshreg

#endif //NEWMESHREG_FEATURESPACE_H
