#include "DiscreteGroupModel.h"

namespace newmeshreg {

void DiscreteGroupModel::initialize_pairs() {

    m_num_pairs = 0;
    for (int n = 0; n < m_num_subjects; n++)
        for (int i = 0; i < control_grid_size; i++)
            for (int n2 = n+1; n2 < m_num_subjects; n2++)
                m_num_pairs++;
    pairs = new int[2 * m_num_pairs];
}

void DiscreteGroupModel::estimate_pairs() {

    int pair = 0;
    std::vector<std::shared_ptr<newresampler::Octree>> cp_grid_trees(m_num_subjects);

    #pragma omp parallel for num_threads(_nthreads)
    for (int subject = 0; subject < m_num_subjects; ++subject)
        cp_grid_trees[subject] = std::make_shared<newresampler::Octree>(m_controlmeshes[subject]);

    for (int subject = 0; subject < m_num_subjects; subject++)
        for (int vertex = 0; vertex < control_grid_size; vertex++) {
            newresampler::Point CP = m_controlmeshes[subject].get_coord(vertex);
            for (int n2 = subject + 1; n2 < m_num_subjects; n2++) {
                int node_ids[2] = {subject * control_grid_size + vertex,
                                   cp_grid_trees[n2]->get_closest_vertex_ID(CP) + n2 * control_grid_size };
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
    constexpr int vertex_per_triangle = 3;

    triplets = new int[vertex_per_triangle * m_num_triplets];

    #pragma omp parallel for num_threads(_nthreads)
    for (int subject = 0; subject < m_num_subjects; subject++)
        for (int triangle = 0; triangle < num_triangles; triangle++)
        {
            int node_ids[vertex_per_triangle] =
                  {m_controlmeshes[subject].get_triangle_vertexID(triangle, 0) + subject * control_grid_size,
                   m_controlmeshes[subject].get_triangle_vertexID(triangle, 1) + subject * control_grid_size,
                   m_controlmeshes[subject].get_triangle_vertexID(triangle, 2) + subject * control_grid_size };
            std::sort(std::begin(node_ids), std::end(node_ids));
            triplets[subject * num_triangles * vertex_per_triangle + triangle * vertex_per_triangle    ] = node_ids[0];
            triplets[subject * num_triangles * vertex_per_triangle + triangle * vertex_per_triangle + 1] = node_ids[1];
            triplets[subject * num_triangles * vertex_per_triangle + triangle * vertex_per_triangle + 2] = node_ids[2];
        }
}

void DiscreteGroupModel::get_rotations(std::vector<NEWMAT::Matrix>& ROT) {

    ROT.clear();
    ROT.resize(m_num_subjects * control_grid_size);
    const newresampler::Point ci = m_samplinggrid.get_coord(m_centroid);

    #pragma omp parallel for num_threads(_nthreads)
    for (int subject = 0; subject < m_num_subjects; subject++)
        for (int vertex = 0; vertex < control_grid_size; vertex++)
            ROT[subject * control_grid_size + vertex] = estimate_rotation_matrix(ci, m_controlmeshes[subject].get_coord(vertex));
}

void DiscreteGroupModel::Initialize(const newresampler::Mesh& controlgrid) {

    control_grid_size = controlgrid.nvertices();
    m_num_nodes = control_grid_size * m_num_subjects;
    m_maxs_dist = 0.5 * controlgrid.calculate_MaxVD();

    m_controlmeshes.clear();
    m_controlmeshes.resize(m_num_subjects, controlgrid);

    datameshtrees.clear();
    datameshtrees.resize(m_num_subjects);

    #pragma omp parallel for num_threads(_nthreads)
    for(int subject = 0; subject < m_num_subjects; ++subject)
        datameshtrees[subject] = std::make_shared<newresampler::Octree>(m_datameshes[subject]);

    std::vector<NEWMAT::ColumnVector> spacings(m_num_subjects);
    #pragma omp parallel for num_threads(_nthreads)
    for(int subject = 0; subject < m_num_subjects; subject++)
    {
        NEWMAT::ColumnVector vMAXmvd(control_grid_size);
        vMAXmvd = 0;
        for (int vertex = 0; vertex < control_grid_size; vertex++)
        {
            newresampler::Point CP = m_controlmeshes[subject].get_coord(vertex);
            for (auto it = m_controlmeshes[subject].nbegin(vertex); it != m_controlmeshes[subject].nend(vertex); it++)
            {
                double dist = 2*RAD * asin((CP - m_controlmeshes[subject].get_coord(*it)).norm() / (2 * RAD));
                if(dist > vMAXmvd(vertex + 1)) vMAXmvd(vertex + 1) = dist;
            }
        }
        spacings[subject] = vMAXmvd;
    }

    //---GET BETWEEN MESH GROUPINGS---//
    initialize_pairs();
    estimate_pairs();

    estimate_triplets();

    initLabeling();
    Initialize_sampling_grid();

    costfct->set_meshes(m_template, m_datameshes[0], controlgrid, m_num_subjects);
    costfct->set_trees(datameshtrees);
    costfct->set_group_spacings(spacings);
    m_iter = 1;
}

void DiscreteGroupModel::setupCostFunction() {

    resetLabeling(); // initialise label array to zero

    for (int n = 0; n < m_num_subjects; n++)
        costfct->reset_CPgrid(m_controlmeshes[n], n);

    costfct->set_iter(m_iter);

    if(m_iter % 2 == 0) m_labels = m_samples;
    else m_labels = m_barycentres;

    m_num_labels = m_labels.size();

    get_rotations(m_ROT);
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
