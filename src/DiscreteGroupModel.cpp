#include "DiscreteGroupModel.h"

namespace newmeshreg {

void DiscreteGroupModel::initialize_pairs() {

    m_num_pairs = 0;
    for (int n = 0; n < m_num_subjects; n++)
        for (int i = 0; i < m_controlmeshes[n].nvertices(); i++)
            for (int n2 = n + 1; n2 < m_num_subjects; n2++)
                //if (n2 > n)
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
                    pairs[2 * pair    ] = node_ids[0];
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

    between_subject_pairs.clear();
    between_subject_pairs.resize(m_num_subjects, std::vector<std::vector<int>>(m_num_subjects,std::vector<int> (control_grid_size,-1)));

    std::vector<std::shared_ptr<newresampler::Octree>> cp_grid_trees(m_num_subjects);

    for (int i = 0; i < m_num_subjects; ++i)
        cp_grid_trees[i] = std::make_shared<newresampler::Octree>(m_controlmeshes[i]);

    for (int n = 0; n < m_num_subjects; ++n)
        for (int i = 0; i < control_grid_size; ++i) {
            newresampler::Point tmp = m_controlmeshes[n].get_coord(i);
            for (int subject = 0; subject < m_num_subjects; ++subject)
                if (n != subject)
                    between_subject_pairs[n][subject][i] = cp_grid_trees[subject]->get_closest_vertex_ID(tmp) * subject;
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
    datameshtrees.resize(m_num_subjects);

    if(m_debug)
    {
        m_template.save("TEMPLATE_res" + std::to_string(m_CPres) + ".surf.gii");
        m_template_LR.save("TEMPLATE_res" + std::to_string(m_CPres) + ".surf.gii");
    }

    initialize_pairs();
    estimate_triplets();

    MVD = controlgrid.calculate_MeanVD();
    m_maxs_dist = 0.5 * controlgrid.calculate_MaxVD();

    initLabeling();
    Initialize_sampling_grid();

    for(int subject = 0; subject < m_num_subjects; ++subject)
        datameshtrees[subject] = std::make_shared<newresampler::Octree>(m_datameshes[subject]);

    costfct->set_meshes(m_template, m_datameshes[0], controlgrid, m_datameshes.size());
    costfct->set_trees(datameshtrees);

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
        std::cout << " initialize cost function " << m_iter << " m_num_triplets " << m_num_triplets << std::endl;

    costfct->setPairs(pairs);
    costfct->setTriplets(triplets);

    costfct->initialize(m_num_nodes, m_num_labels, m_num_pairs, m_num_triplets);
    costfct->get_source_data();

    if(m_verbosity)
        std::cout << " numpoints " << m_num_nodes << " m_num_labels " << m_num_labels << " m_num_pairs " << m_num_pairs << std::endl;

    m_iter++;
}

} //namespace newmeshreg
