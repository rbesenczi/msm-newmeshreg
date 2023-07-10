#include "DiscreteGroupCostFunction.h"

namespace newmeshreg {

void DiscreteGroupCostFunction::set_parameters(myparam& p) {
    myparam::iterator it;
    it=p.find("lambda");_reglambda=boost::get<float>(it->second);
    it=p.find("lambda_pairs");_lambdapairs=boost::get<float>(it->second);
    it=p.find("set_lambda_pairs");_setpairs=boost::get<bool>(it->second);
    it=p.find("simmeasure");_simmeasure=boost::get<int>(it->second); sim.set_simval(_simmeasure);
    it=p.find("verbosity");_verbosity=boost::get<bool>(it->second);
    it=p.find("shearmodulus");_mu=boost::get<float>(it->second);
    it=p.find("bulkmodulus");_kappa=boost::get<float>(it->second);
    it=p.find("sigma_in");_sigma=boost::get<float>(it->second);
    //it=p.find("quartet");_quadcost=boost::get<bool>(it->second);
}

void DiscreteGroupCostFunction::get_spacings() {

    NEWMAT::ColumnVector vMAXmvd(_CONTROLMESHES[0].nvertices());
    SPACINGS.clear();

    for(int n = 0; n < num_subjects; n++)
    {
        vMAXmvd = 0;
        for (int k = 0; k < _CONTROLMESHES[n].nvertices(); k++)
        {
            newresampler::Point CP = _CONTROLMESHES[n].get_coord(k);
            for (auto it = _CONTROLMESHES[n].nbegin(k); it != _CONTROLMESHES[n].nend(k); it++)
            {
                double dist = 2 * RAD * asin((CP -_CONTROLMESHES[n].get_coord(*it)).norm() / (2*RAD));
                if(dist > vMAXmvd(k+1)) vMAXmvd(k+1) = dist;
            }
        }
        SPACINGS.push_back(vMAXmvd);
    }
}

void DiscreteGroupCostFunction::initialize(int numNodes, int numLabels, int numPairs, int numTriplets,
                                           int numQuartets) {

    DiscreteCostFunction::initialize(numNodes,numLabels,numPairs,numTriplets,numQuartets);

    define_template_patches();
    resample_to_template();
    get_spacings();
    resample_patches();
/*
    if (_quadcost && !_setpairs)
        _lambdapairs = (double) m_num_quartets / (double) m_num_pairs;*/
}

void DiscreteGroupCostFunction::define_template_patches() {

    std::vector<std::vector<int>> TOTALCELLS(m_num_nodes,std::vector<int>());
    std::vector<std::vector<int>> TEMPLATE_CELLS;
    TEMPLATEPTS.clear();

    for (int n = 0; n < m_num_nodes; n++)
    {
        int mesh_ID = std::floor(n / VERTICES_PER_SUBJ);
        int ind = n - mesh_ID * VERTICES_PER_SUBJ;
        newresampler::Point p = _CONTROLMESHES[mesh_ID].get_coord(ind);
        TEMPLATE_CELLS.push_back(_targetrel->return_cell_group(p,1.75*asin(MVD_LR/RAD)));
        TEMPLATEPTS.emplace_back(std::vector<int>());
    }

    // make patch constant across all node pairs
    for (int pair = 0; pair < m_num_pairs; pair++)
        for (int i = 0; i < 2; i++)
            for (int j = 0; j < TEMPLATE_CELLS[_pairs[2 * pair + i]].size(); j++)
                TOTALCELLS[_pairs[2 * pair + i]].push_back(TEMPLATE_CELLS[_pairs[2 * pair + i]][j]);

    for (int n = 0; n < m_num_nodes; n++)
    {
        std::sort(TOTALCELLS[n].begin(),TOTALCELLS[n].end());
        TOTALCELLS[n].erase( unique(TOTALCELLS[n].begin(),TOTALCELLS[n].end()),TOTALCELLS[n].end());

        for(int i = 0; i < TOTALCELLS[n].size(); i++)
        {
            std::vector<int> PTS = _targetrel->get_cell_members(TOTALCELLS[n][i]);
            TEMPLATEPTS[n].insert(TEMPLATEPTS[n].end(),PTS.begin(),PTS.end());
        }
    }
}

void DiscreteGroupCostFunction::resample_to_template() {

    RESAMPLEDDATA.clear();
    for(int n = 0; n < num_subjects; n++)
    {
        RESAMPLEDDATA.push_back(FEAT->get_data_matrix(n));
        _TEMPLATE.set_pvalues(RESAMPLEDDATA[n]);
        newresampler::metric_resample(_DATAMESHES[n],_TEMPLATE);
        //R.resampledata(_DATAMESHES[n],_TEMPLATE,RESAMPLEDDATA[n],0.0,_targetrel);
    }
}
/*
double DiscreteGroupCostFunction::computeQuartetCost(int quartet, int labelA, int labelB, int labelC, int labelD) {

    double cost=0;
    std::vector<int> indices(4,0);
    std::vector<int> rows;
    int mesh_ID, index = 0;
    std::map<int,int> testsubjects;
    bool found;

    indices[0]=labelA+_quartets[4*quartet]*_labels.size();
    indices[1]=labelB+_quartets[4*quartet+1]*_labels.size();
    indices[2]=labelC+_quartets[4*quartet+2]*_labels.size();
    indices[3]=labelD+_quartets[4*quartet+3]*_labels.size();

    for(int i = 0; i < 4; i++)
    {
        mesh_ID = std::floor(_quartets[4*quartet+i]/VERTICES_PER_SUBJ);
        testsubjects[mesh_ID] = indices[i];
    }

    for (auto iter = PATCHDATA[indices[0]].begin(); iter != PATCHDATA[indices[0]].end(); ++iter)
    {
        found = true;
        for (int i = 1; i < 4; i++)
            if (PATCHDATA[indices[i]].find(iter->first) == PATCHDATA[indices[i]].end())
                found = false;
        if (found)
            rows.push_back(iter->first);
    }

    NEWMAT::Matrix Groupdata(rows.size(),num_subjects);
    Groupdata = 0;

    if(rows.size() > 0)
    {
        for (int n = 0; n < num_subjects; n++)
            for (int r = 0; r < rows.size(); r++)
                // first add all the test subjects then add remaining subjects?
                if (testsubjects.find(n) != testsubjects.end())
                    Groupdata(r + 1, n + 1) = PATCHDATA[testsubjects[n]].find(rows[r])->second;
                else
                    Groupdata(r + 1, n + 1) = RESAMPLEDDATA[n](1, rows[r] + 1);


        NEWMAT::DiagonalMatrix eigenvals;
        SVD(Groupdata,eigenvals);

        for (int i = 1; i <= eigenvals.Nrows(); i++)
            cost += eigenvals(i);
    }
    else
        cost = 1e6;

    return _lambdapairs * cost;
}
*/
double DiscreteGroupCostFunction::computeTripletCost(int triplet, int labelA, int labelB, int labelC) {

    int meshID = floor(triplet/TRIPLETS_PER_SUBJ);
    int m0 = _triplets[3*triplet]-meshID*VERTICES_PER_SUBJ;
    int m1 = _triplets[3*triplet+1]-meshID*VERTICES_PER_SUBJ;
    int m2 = _triplets[3*triplet+2]-meshID*VERTICES_PER_SUBJ;

    newresampler::Triangle TRI((*ROTATIONS)[_triplets[3*triplet  ]] * _labels[labelA],
                               (*ROTATIONS)[_triplets[3*triplet+1]] * _labels[labelB],
                               (*ROTATIONS)[_triplets[3*triplet+2]] * _labels[labelC],0);
    newresampler::Triangle TRI_ORIG(_ORIG.get_coord(m0),
                                    _ORIG.get_coord(m1),
                                    _ORIG.get_coord(m2),0);
    newresampler::Triangle TRI_noDEF(_CONTROLMESHES[meshID].get_coord(m0),
                                     _CONTROLMESHES[meshID].get_coord(m1),
                                     _CONTROLMESHES[meshID].get_coord(m2),0);

    if((TRI.normal() | TRI_noDEF.normal()) < 0) { return 1e7; }

    return _reglambda * MISCMATHS::pow(calculate_triangular_strain(TRI_ORIG,TRI,_mu,_kappa),_rexp);
}

double DiscreteGroupCostFunction::computePairwiseCost(int pair, int labelA, int labelB) {
    return sim.corr(PATCHDATA[labelA+_pairs[2*pair]*_labels.size()], PATCHDATA[labelB+_pairs[2*pair+1]*_labels.size()]);
}

void DiscreteGroupCostFunction::resample_patches() {

    std::vector<bool> withinCPrange(m_num_nodes*_labels.size(),false);
    PATCHDATA.clear();
    PATCHDATA.resize(m_num_nodes*_labels.size(), std::map<int,float>());

    for (int n = 0; n < m_num_nodes; n++)
    {
        int mesh_ID = floor(n/VERTICES_PER_SUBJ);
        int node = n - mesh_ID * VERTICES_PER_SUBJ;
        for (int lab = 0; lab < _labels.size(); lab++)
            if ((_labels[lab] - _labels[0]).norm() <= 0.5 * SPACINGS[mesh_ID](node + 1))
                withinCPrange[lab + n * _labels.size()] = true;

    }
    resampler_worker_function(0,m_num_nodes*_labels.size(),withinCPrange);
}

void DiscreteGroupCostFunction::resampler_worker_function(int begin, int end, const std::vector<bool>& INrange) {

    for (int n = begin; n < end; n++)
    {
        int num = floor(n/_labels.size());
        int mesh_ID = floor(num/VERTICES_PER_SUBJ);
        int node = num - mesh_ID * VERTICES_PER_SUBJ;
        int labnum = n - num * _labels.size();

        newresampler::Point newpt = (*ROTATIONS)[num] * _labels[labnum];

        if (INrange[labnum + num * _labels.size()])
            PATCHDATA[labnum + num * _labels.size()] = resample_onto_template(node, mesh_ID, newpt, TEMPLATEPTS[num]);
    }
}

std::map<int,float> DiscreteGroupCostFunction::resample_onto_template(int node, int ID, const newresampler::Point& newCP,
                                                                      const std::vector<int>& PTS) {
    std::map<int,float> sampledata;
    std::map<int,float> weights;

    for(const int& i : PTS)
    {
        sampledata[i] = 0.0;
        weights[i] = 0.0;
    }

    for (auto j = _CONTROLMESHES[ID].tIDbegin(node); j != _CONTROLMESHES[ID].tIDend(node); j++)
    {
        newresampler::Point newCP0,newCP1,newCP2;
        newresampler::Triangle tri = _CONTROLMESHES[ID].get_triangle(*j);
        newresampler::Point CP0 = tri.get_vertex_coord(0);
        newresampler::Point CP1 = tri.get_vertex_coord(1);
        newresampler::Point CP2 = tri.get_vertex_coord(2);
        int id0 = tri.get_vertex_no(0);
        int id1 = tri.get_vertex_no(1);
        int id2 = tri.get_vertex_no(2);

        if(id0==node)
        {
            newCP0 = newCP;
            newCP1 = CP1;
            newCP2 = CP2;
        }
        else if(id1==node)
        {
            newCP0 = CP0;
            newCP1 = newCP;
            newCP2 = CP2;
        }
        else if(id2==node)
        {
            newCP0 = CP0;
            newCP1 = CP1;
            newCP2 = newCP;
        }
        else
        {
            std::cout << node << " triangle node not found " << " id0" << id0 << " " << id1 << " " << id2 << std::endl;
            exit(1);
        }

        for(int i = 1; i <= _sourcerel.Nrows(*j+1); i++)
        {
            int sourceind = _sourcerel(i,*j+1);
            // move source using barycentric interpolation
            newresampler::Point SP = _DATAMESHES[ID].get_coord(sourceind-1);
            newresampler::project_point(CP0, CP1, CP2,SP);
            newresampler::Point tmp = barycentric(CP0,CP1,CP2,SP,newCP0,newCP1,newCP2);
            tmp.normalize(); //ASSUMING SPHERE"S ORIGIN IS 0,0,0 (which  we can project from the face to the surface of the sphere this way
            tmp = tmp * RAD;

            // resample onto template using gaussian
            std::vector<std::pair<float,int>> NEIGHBOURS = _targetrel->update_RELATIONS_for_ind(tmp);
            for (const auto& s : NEIGHBOURS)
            {
                double dist = (tmp - _TEMPLATE.get_coord(s.second)).norm();
                double g_weight = exp(-(dist*dist)/(2*_sigma*_sigma));
                sampledata[s.second] += g_weight * FEAT->get_data_val(1,sourceind,ID);
                weights[s.second] += g_weight;
            }
        }
    }

    for(int i : PTS)
    {
        if(weights[i]==0.0)
            sampledata[i]=RESAMPLEDDATA[ID](1,i+1);
        else
            sampledata[i] /= weights[i];
    }

    return sampledata;
}

} //namespace newmeshreg
