#include "DiscreteGroupModel.h"

namespace newmeshreg {

void DiscreteGroupModel::initialize_pairs() {

    m_num_pairs = 0;
    for (int n = 0; n < m_num_subjects; n++)
        for (int i = 0; i < m_controlmeshes[n].nvertices(); i++)
            for (int n2 = n+1; n2 < m_num_subjects; n2++)
                m_num_pairs++;
    pairs = new int[2 * m_num_pairs];
}

void DiscreteGroupModel::estimate_pairs() {

    int pair = 0;

    std::vector<std::shared_ptr<newresampler::Octree>> cp_grid_trees(m_num_subjects);

    #pragma omp parallel for num_threads(_nthreads)
    for (int i = 0; i < m_num_subjects; ++i)
        cp_grid_trees[i] = std::make_shared<newresampler::Octree>(m_controlmeshes[i]);

    for (int n = 0; n < m_num_subjects; n++)
        for (int i = 0; i < control_grid_size; i++) {
            newresampler::Point tmp = m_controlmeshes[n].get_coord(i);
            for (int n2 = n + 1; n2 < m_num_subjects; n2++)
            {
                int node_ids[2] = {n * control_grid_size + i,
                                   cp_grid_trees[n2]->get_closest_vertex_ID(tmp) + n2 * control_grid_size };
                std::sort(std::begin(node_ids), std::end(node_ids));
                pairs[2*pair  ] = node_ids[0];
                pairs[2*pair+1] = node_ids[1];
                pair++;
            }
        }
}

void DiscreteGroupModel::estimate_triplets() {

    int num_triangles = m_controlmeshes[0].ntriangles();
    m_num_triplets = m_num_subjects * num_triangles;

    triplets = new int[3 * m_num_triplets];

    #pragma omp parallel for num_threads(_nthreads)
    for (int n = 0; n < m_num_subjects; n++)
        for (int i = 0; i < num_triangles; i++)
        {
            int node_ids[3] = { m_controlmeshes[n].get_triangle_vertexID(i, 0) + n * control_grid_size,
                                m_controlmeshes[n].get_triangle_vertexID(i, 1) + n * control_grid_size,
                                m_controlmeshes[n].get_triangle_vertexID(i, 2) + n * control_grid_size };
            std::sort(std::begin(node_ids), std::end(node_ids));
            triplets[3 * n * num_triangles + 3 * i    ] = node_ids[0];
            triplets[3 * n * num_triangles + 3 * i + 1] = node_ids[1];
            triplets[3 * n * num_triangles + 3 * i + 2] = node_ids[2];
        }
}

void DiscreteGroupModel::get_rotations(std::vector<NEWMAT::Matrix>& ROT) {

    ROT.clear();
    ROT.resize(m_num_subjects * control_grid_size);
    const newresampler::Point& ci = m_samplinggrid.get_coord(m_centroid);

    #pragma omp parallel for num_threads(_nthreads)
    for (int n = 0; n < m_num_subjects; n++)
        for (int k = 0; k < m_controlmeshes[n].nvertices(); k++)
            ROT[n * control_grid_size + k] = estimate_rotation_matrix(ci, m_controlmeshes[n].get_coord(k));
}

void DiscreteGroupModel::Initialize(const newresampler::Mesh& controlgrid) {

    control_grid_size = controlgrid.nvertices();
    m_num_nodes = control_grid_size * m_num_subjects;
    m_controlmeshes.clear();
    m_controlmeshes.resize(m_num_subjects, controlgrid);

    initialize_pairs();
    estimate_triplets();

    m_maxs_dist = 0.5 * controlgrid.calculate_MaxVD();

    initLabeling();
    Initialize_sampling_grid();

    costfct->set_meshes(m_template, m_datameshes[0], controlgrid, m_num_subjects);
    m_iter = 1;
}

void DiscreteGroupModel::setupCostFunction() {

    resetLabeling(); // initialise label array to zero

    for (int n = 0; n < m_num_subjects; n++)
        costfct->reset_CPgrid(m_controlmeshes[n], n);

    costfct->set_iter(m_iter);

    datameshtrees.clear();
    datameshtrees.resize(m_num_subjects);

    #pragma omp parallel for num_threads(_nthreads)
    for(int subject = 0; subject < m_num_subjects; ++subject)
        datameshtrees[subject] = std::make_shared<newresampler::Octree>(m_datameshes[subject]);
    costfct->set_trees(datameshtrees);

    //---GET BETWEEN MESH GROUPINGS---//
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

    if(m_verbosity)
        std::cout << " numpoints " << m_num_nodes << " m_num_labels " << m_num_labels << " m_num_pairs " << m_num_pairs << std::endl;

    m_iter++;
}

} //namespace newmeshreg
