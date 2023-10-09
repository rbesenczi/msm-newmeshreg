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

void DiscreteGroupModel::get_rotated_meshes() {
    rotated_datameshes.clear();
    resampled_data.clear();
    rotated_datameshes.resize(m_num_subjects*m_num_labels);
    resampled_data.resize(m_num_subjects*m_num_labels);
    for(int subject = 0; subject < m_num_subjects; subject++) {
        newresampler::Mesh rotated_mesh = m_datameshes[subject];
        rotated_mesh.set_pvalues(FEAT->get_data_matrix(subject));
        #pragma omp parallel for num_threads(_nthreads)
        for(int label = 0; label < m_num_labels; label++) {
            for (int vertex = 0; vertex < control_grid_size; vertex++) {
                const newresampler::Point CP = m_controlmeshes[subject].get_coord(vertex);
                for (int datapoint = 0; datapoint < m_datameshes[subject].nvertices(); datapoint++) {
                    const newresampler::Point SP = m_datameshes[subject].get_coord(datapoint);
                    if (((2 * RAD * asin((CP - SP).norm() / (2 * RAD))) < range * spacings[subject](vertex + 1)))
                        rotated_mesh.set_coord(datapoint,
                                               SP*estimate_rotation_matrix(CP, m_ROT[subject*control_grid_size+vertex]*m_labels[label]));
                }
            }
            newresampler::Mesh tmp = newresampler::metric_resample(rotated_mesh, m_template);
            rotated_datameshes.at(subject*m_num_labels+label) = tmp;
            resampled_data.at(subject*m_num_labels+label) = tmp.get_pvalues();
        }
    }
}

void DiscreteGroupModel::get_patch_data() {
    patch_data.clear();
    patch_data.resize(control_grid_size * m_num_subjects * m_num_labels);
    for (int subject = 0; subject < m_num_subjects; subject++)
        #pragma omp parallel for num_threads(_nthreads)
        for (int vertex = 0; vertex < control_grid_size; vertex++) {
            for (int label = 0; label < m_num_labels; label++)  {
                std::map<int, double> patchdata;
                const newresampler::Point CP = m_ROT[subject * control_grid_size + vertex] * m_labels[label];
                for (int datapoint = 0; datapoint < m_template.nvertices(); datapoint++)
                    if (((2 * RAD * asin((CP - m_template.get_coord(datapoint)).norm() / (2 * RAD))) < range * spacings[subject](vertex + 1)))
                        patchdata[datapoint] = resampled_data[subject*m_num_labels+label](1,datapoint + 1);
                patch_data[subject * control_grid_size * m_num_labels + vertex * m_num_labels + label] = patchdata;
            }
        }
}

void DiscreteGroupModel::Initialize(const newresampler::Mesh& controlgrid) {

    control_grid_size = controlgrid.nvertices();
    m_num_nodes = control_grid_size * m_num_subjects;
    m_maxs_dist = 0.5 * controlgrid.calculate_MaxVD();

    m_controlmeshes.clear();
    m_controlmeshes.resize(m_num_subjects, controlgrid);

    targettree = std::make_shared<newresampler::Octree>(m_template);

    spacings.clear();
    spacings.resize(m_num_subjects);
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

    m_labels = m_samples;
    m_num_labels = m_labels.size();

    get_rotations(m_ROT);
    //get_rotated_meshes();
    //get_patch_data();

    costfct->set_labels(m_labels,m_ROT);
    costfct->set_meshes(m_template, m_datameshes[0], controlgrid, m_num_subjects);
    costfct->set_targettree(targettree);
    costfct->set_group_spacings(spacings);
    //costfct->set_rotated_meshes(rotated_datameshes, resampled_data);
    //costfct->set_patch_data(patch_data);
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
    get_rotated_meshes();
    get_patch_data();

    //---INIT---//
    if(m_verbosity)
        std::cout << " initialize cost function " << m_iter << " m_num_triplets " << m_num_triplets << std::endl;

    costfct->setPairs(pairs);
    costfct->setTriplets(triplets);
    costfct->set_patch_data(patch_data);

    costfct->initialize(m_num_nodes, m_num_labels, m_num_pairs, m_num_triplets);

    if(m_verbosity)
        std::cout << " numpoints " << m_num_nodes << " m_num_labels " << m_num_labels << " m_num_pairs " << m_num_pairs << std::endl;

    m_iter++;
}

} //namespace newmeshreg
