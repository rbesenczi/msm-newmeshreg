#include "DiscreteGroupModel.h"

namespace newmeshreg {

void DiscreteGroupModel::initialize_pairs() {

    m_num_pairs = 0;
    for (int n = 0; n < m_num_subjects; n++)
        for (int i = 0; i < m_controlmeshes[n].nvertices(); i++)
            for (int n2 = 0; n2 < m_num_subjects; n2++)
                if (n2 > n)
                    m_num_pairs++;
}

void DiscreteGroupModel::estimate_pairs() {

    int pair = 0;
    pairs = new int[m_num_pairs*2];

    for (int n = 0; n < m_num_subjects; n++)
        for (int i = 0; i < m_controlmeshes[n].nvertices(); i++)
            for (int n2 = 0; n2 < m_num_subjects; n2++)
                if (n2 > n)
                {
                    int node_ids[2] = { n * control_grid_size + i,between_subject_pairs[n][n2][i] };
                    std::sort(std::begin(node_ids), std::end(node_ids));
                    pairs[2 * pair] = node_ids[0];
                    pairs[2 * pair + 1] = node_ids[1];
                    pair++;
                }
}

void DiscreteGroupModel::estimate_triplets() {

    m_num_triplets = m_num_subjects * m_controlmeshes[0].ntriangles();

    triplets = new int[m_num_triplets*3];

    for (int n = 0; n < m_num_subjects; n++)
        for (int i = 0; i < m_controlmeshes[n].ntriangles(); i++)
        {
            int node_ids[3] = { m_controlmeshes[n].get_triangle_vertexID(i, 0) + n * m_controlmeshes[n].nvertices(),
                                m_controlmeshes[n].get_triangle_vertexID(i, 1) + n * m_controlmeshes[n].nvertices(),
                                m_controlmeshes[n].get_triangle_vertexID(i, 2) + n * m_controlmeshes[n].nvertices() };
            std::sort(std::begin(node_ids), std::end(node_ids));
            triplets[3 * i + n * 3 * m_controlmeshes[n].ntriangles()    ] = node_ids[0];
            triplets[3 * i + n * 3 * m_controlmeshes[n].ntriangles() + 1] = node_ids[1];
            triplets[3 * i + n * 3 * m_controlmeshes[n].ntriangles() + 2] = node_ids[2];
        }
}

void DiscreteGroupModel::get_between_subject_pairs() {

    //TODO need to work on m_TEMPLATE_LR_ALL_RELATIONS...

    double ang = 2 * asin(MVD / RAD);
    between_subject_pairs.clear();
    between_subject_pairs.resize(m_num_subjects, std::vector<std::vector<int>>(m_num_subjects,std::vector<int> (control_grid_size,-1)));

    // create one giant mesh for all subject control grids
    std::vector<std::shared_ptr<newresampler::Mpoint>> ALLPOINTS = m_controlmeshes[0].get_points();
    for (int n = 1; n < m_num_subjects; n++)
    {
        std::vector<std::shared_ptr<newresampler::Mpoint>> tmppoints = m_controlmeshes[n].get_points();
        ALLPOINTS.insert(ALLPOINTS.end(),tmppoints.begin(),tmppoints.end());
    }

    m_TEMPLATE_LR_ALL_RELATIONS = std::shared_ptr<RELATIONS>(new RELATIONS(m_template_LR, ALLPOINTS, ang));

    // for all subjects and all vertices find the closest between mesh neighbours
    for(int n = 0; n < m_num_subjects; n++)
    {
        double angtmp = ang;
        int i = 1;
        while (i <= control_grid_size)
        {
            int found = 0;
            m_TEMPLATE_LR_ALL_RELATIONS->update_RELATIONS_for_ind(i, m_controlmeshes[n], angtmp);
            for (int r = 1; r <= m_TEMPLATE_LR_ALL_RELATIONS->Nrows(i); r++)
            {
                int ind = (*m_TEMPLATE_LR_ALL_RELATIONS)(r,i)-1;
                int mesh_ID = floor(ind / control_grid_size);
                if(mesh_ID != n && between_subject_pairs[n][mesh_ID][i-1] == -1)
                {
                    between_subject_pairs[n][mesh_ID][i-1] = ind;
                    found++;
                    if(found == m_num_subjects - 1) break;
                }
            }
            if (found != m_num_subjects - 1) {
                angtmp = 2 * ang; // check that neighbours have been found between all other meshes
                break;
            }
            else
                i++;
        }
    }
}

void DiscreteGroupModel::get_rotations(std::vector<NEWMAT::Matrix>& ROT) {

    ROT.clear();
    const newresampler::Point& ci = m_samplinggrid.get_coord(m_centroid);

    for (int n = 0; n < m_num_subjects; n++)
        for (int k = 0; k < m_controlmeshes[n].nvertices(); k++)
            ROT.emplace_back(estimate_rotation_matrix(ci, m_controlmeshes[n].get_coord(k)));
}

void DiscreteGroupModel::Initialize(const newresampler::Mesh& controlgrid) {

    control_grid_size = controlgrid.nvertices();
    m_num_nodes = control_grid_size * m_num_subjects;
    m_controlmeshes.clear();
    m_controlmeshes.resize(m_num_subjects, controlgrid);
    m_template_LR = controlgrid;

    if(m_debug)
    {
        m_template.save("TEMPLATE_res" + std::to_string(m_CPres) + ".surf.gii");
        m_template_LR.save("TEMPLATE_res" + std::to_string(m_CPres) + ".surf.gii");
    }

    /*
    for(int n = 0; n < m_num_subjects; n++)
        m_controlmeshes.push_back(controlgrid);
    */
    //---INITIALIZE REGULARISATION TRIPLETS---//
    initialize_pairs();
    estimate_triplets();

    //---CALCULATE FIXED GRID SPACINGS---//
    MVD = controlgrid.calculate_MeanVD();
    m_maxs_dist = 0.4 * controlgrid.calculate_MaxVD();

    //---INITIALIAZE LABEL GRID---//
    initLabeling();
    Initialize_sampling_grid();

    //---INITIALIZE NEIGHBOURHOODS---//
    //m_cp_neighbourhood=std::shared_ptr<RELATIONS>(new RELATIONS(m_DATAMESHES[0], controlgrid, 2 * asin(MVD / RAD)));
    // as source mesh moves with control grid the relationships are constant for all meshes
    //m_inputrel=std::shared_ptr<RELATIONS>(new RELATIONS(m_DATAMESHES[0],m_TEMPLATE,2*asin(MVD/RAD)));
    //m_inputtree = std::make_shared<newresampler::Octree>(m_template);
    //m_cp_neighbourhood->update_RELATIONS(m_DATAMESHES[0]);
    //costfct->set_relations(m_cp_neighbourhood,m_inputrel);
    //costfct->set_octrees(m_inputtree);

    costfct->set_meshes(m_template, m_datameshes[0], controlgrid, m_datameshes.size());

    m_iter = 1;
}

void DiscreteGroupModel::setupCostFunction() {

    resetLabeling(); // initialise label array to zero

    //---use geodesic distances---//
    for (int n = 0; n < m_num_subjects; n++)
        costfct->reset_CPgrid(m_controlmeshes[n], n);

    costfct->set_iter(m_iter);

    //---GET BETWEEN MESH GROUPINGS---//
    get_between_subject_pairs();

    estimate_pairs();

    //---GET LABEL SPACE---//
    get_rotations(m_ROT);
    if(m_iter % 2 == 0) m_labels = m_samples;
    else m_labels = m_barycentres;

    m_num_labels = m_labels.size();

    costfct->set_labels(m_labels,m_ROT);

    //---INIT---//
    if(m_verbosity)
        std::cout << " initialize cost function 2" << m_iter << " m_num_triplets " << m_num_triplets << std::endl;

    costfct->setPairs(pairs);
    costfct->setTriplets(triplets);

    costfct->initialize(m_num_nodes, m_num_labels, m_num_pairs, m_num_triplets);

    if(m_verbosity)
        std::cout << " numpoints " << m_num_nodes << " m_num_labels " << m_num_labels << " m_num_pairs " << m_num_pairs << std::endl;

    m_iter++;
}

void DiscreteGroupModel::applyLabeling(int* dlabels) {
    for (int n = 0; n < m_num_subjects; n++)
        for (int i = 0; i < m_controlmeshes[n].nvertices(); i++)
            m_controlmeshes[n].set_coord(i, m_ROT[i + n * m_controlmeshes[n].nvertices()] *
                                            m_labels[dlabels[i + n * m_controlmeshes[n].nvertices()]]);
}

} //namespace newmeshreg
