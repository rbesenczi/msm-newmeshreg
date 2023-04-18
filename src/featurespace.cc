#include "featurespace.h"

using namespace std;

namespace newmeshreg {

featurespace::featurespace(const std::string& datain, const std::string& dataref) {
    _sigma_in.resize(2,5.0);
    _fthreshold.resize(2,0.0);
    CMfile_in.push_back(datain);
    CMfile_in.push_back(dataref);
}

featurespace::featurespace(const std::string& datain, const std::vector<std::string>& datareflist) {
    _sigma_in.push_back(5);
    _fthreshold.resize(2,0.0);
    CMfile_in.push_back(datain);

    for(const auto& i : datareflist)
    {
        CMfile_in.push_back(i);
        _sigma_in.push_back(5);
    }
}

featurespace::featurespace(const std::vector<std::string>& datalist) {
    _fthreshold.resize(2,0.0);

    for(const auto& i : datalist)
    {
        CMfile_in.push_back(i);
        _sigma_in.push_back(5);
    }
}

void featurespace::set_smoothing_parameters(const std::vector<double>& s) {

    _sigma_in.clear();
    if (s.size() != CMfile_in.size())
        if (s.size() == 1)
            for (int i = 0; i < (int) CMfile_in.size(); i++)
                _sigma_in.push_back(s[0]);
        else
            throw newresampler::MeshException("Mewmesh::featurespace:: smoothing sigma size incompatible with data dimensions");
    else
        _sigma_in = s;
}

newresampler::Mesh featurespace::initialize(int ico, vector<newresampler::Mesh> &IN, bool exclude){

    newresampler::Mesh icotmp;

    if (IN.size() != CMfile_in.size())
        throw MeshregException("featurespace::Initialize do not have the same number of datasets and surface meshes");
    else
        DATA.resize(IN.size(), std::shared_ptr<MISCMATHS::BFMatrix>());

    if(ico > 0)
    {
        icotmp = newresampler::make_mesh_from_icosa(ico);
        newresampler::recentre(icotmp);
        newresampler::true_rescale(icotmp, RAD);
    }

    for (unsigned int i = 0; i < IN.size(); i++)
    {
        set_data(CMfile_in[i],DATA[i],IN[i],_issparse);

        if (ico == 0) icotmp = IN[i];

        if (exclude || _cut)
        {
            newresampler::Mesh excl_tmp = newresampler::create_exclusion(IN[i], _fthreshold[0], _fthreshold[1]);
            EXCL.push_back(std::make_shared<newresampler::Mesh>(excl_tmp));
        }
        else
            EXCL.push_back(std::shared_ptr<newresampler::Mesh>());

        newresampler::Mesh tmp = newresampler::metric_resample(IN[i], icotmp, _nthreads, EXCL[i]);

        if (_sigma_in[i] > 0)
            tmp = newresampler::smooth_data(tmp, tmp, _sigma_in[i], _nthreads, EXCL[i]);

        DATA[i] = std::make_shared<MISCMATHS::FullBFMatrix>(tmp.get_pvalues());
    }

    // intensity normalise using histogram matching
    if (_intensitynorm)
        for (unsigned int i = 1; i < IN.size(); i++)
            multivariate_histogram_normalization(*DATA[i], *DATA[0], EXCL[i], EXCL[0], _scale);
            // match input data feature distributions to equivalent in ref, rescale all to first feature in reference if _scale is

    if (_logtransform)
        for (unsigned int i = 0; i < IN.size(); i++)
            log10_transform_and_normalise(*DATA[i]);

    if (_varnorm)
        for (unsigned int i = 0; i < IN.size(); i++)
            variance_normalise(DATA[i], EXCL[i]);

    return icotmp;
}

void featurespace::variance_normalise(std::shared_ptr<MISCMATHS::BFMatrix> &DATA, std::shared_ptr<newresampler::Mesh> &EXCL) {

    vector<vector<double>> _data;

    for (unsigned int i = 1; i <= DATA->Ncols(); i++)
        if (!EXCL || EXCL->get_pvalue(i - 1) > 0)
        {
            vector<double> tmp;
            for (unsigned int k = 1; k <= DATA->Nrows(); k++)
                tmp.push_back(DATA->Peek(k, i));
            _data.push_back(tmp);
        }

    vector<double> mean(_data[0].size(), 0.0),
                   var(_data[0].size(), 0.0);

    for(unsigned int j = 0; j < _data[0].size(); j++)
    {
        for (unsigned int i = 0; i < _data.size(); i++)
        {
            double delta = _data[i][j] - mean[j];
            mean[j] += delta / (i + 1);
            var[j] += delta * (_data[i][j] - mean[j]);
        }

        var[j] = var[j] / (_data.size() - 1);

        for(auto& i : _data)
        {
            i[j] -= mean[j];
            if(var[j] > 0)
                i[j] = i[j] / sqrt(var[j]);
        }
    }

    int ind = 0;
    for (unsigned int i = 1; i <= DATA->Ncols(); i++)
        if (!EXCL || EXCL->get_pvalue(i - 1) > 0)
        {
            for (unsigned int k = 1; k <= DATA->Nrows(); k++)
                DATA->Set(k, i, _data[ind][k - 1]);
            ind++;
        }
}

void featurespace::log10_transform_and_normalise(MISCMATHS::BFMatrix& data){

    for (int i = 1; i <= (int) data.Ncols(); i++)
        for (auto it = data.begin(i); it != data.end(i); ++it)
            data.Set(it.Row(), i, log10(*it + 1));

    for(int i = 1; i <= (int) data.Ncols(); i++)
    {
        double size = 0.0, mean = 0.0;
        for (auto it = data.begin(i); it != data.end(i); ++it)
            mean += (*it);

        mean /= data.Nrows();
        for (auto it = data.begin(i); it != data.end(i); ++it)
            size += (*it - mean) * (*it - mean);

        size = sqrt(size);
        for (auto it = data.begin(i); it != data.end(i); ++it)
            if (*it > 0) data.Set(it.Row(), i, *it / size);
    }
}

} //namespace newmeshreg
