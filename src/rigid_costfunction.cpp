#include "rigid_costfunction.h"

#include <memory>

using namespace std;

namespace newmeshreg {

Rigid_cost_function::Rigid_cost_function(
        const newresampler::Mesh& target,
        const newresampler::Mesh& source,
        std::shared_ptr<featurespace>& features)
        :TARGET(target), SOURCE(source), FEAT(features), simmeasure(3), iters(20), verbosity(false)
{
    // if no cf weighting then set weighting to 1 by default
    inweight.ReSize(1, SOURCE.nvertices());
    refweight.ReSize(1, TARGET.nvertices());
    inweight = 1;
    refweight = 1;
}

void Rigid_cost_function::Initialize(){

    current_sim.ReSize(SOURCE.nvertices());
    MVD = SOURCE.calculate_MeanVD();
    min_sigma = MVD;
    sim.Resize(TARGET.nvertices(), SOURCE.nvertices());  // creates sparse sim kernel to fill with similarities of nearest neighbours only
    sim.set_reference(FEAT->get_reference_data());
    sim.set_input(FEAT->get_input_data());
    sim.Initialize(simmeasure);
    nbh = std::make_shared<Neighbourhood>();
    nbh->update(SOURCE, TARGET, 2 * asin(4 * MVD / (2 * RAD)));
    sim.set_neighbourhood(nbh);
    targettree = std::make_shared<newresampler::Octree>(TARGET);
    update_similarity();
}

void Rigid_cost_function::set_parameters(myparam& PAR){
    myparam::iterator it;
    it = PAR.find("iters"); iters = boost::get<int>(it->second);
    it = PAR.find("simmeasure"); simmeasure = boost::get<int>(it->second);
    it = PAR.find("verbosity"); verbosity = boost::get<bool>(it->second);
    it = PAR.find("stepsize"); stepsize = boost::get<float>(it->second);
    it = PAR.find("gradsampling"); spacing = boost::get<float>(it->second);
}

NEWMAT::ColumnVector Rigid_cost_function::WLS_simgradient(const newresampler::Tangs& tangs, int index, const std::vector<int>& querypoints) {
// weighted least squares similarity gradient estimation

    double x11, x21, y11, y21, SUM = 0.0, JPsim = 0.0;
    NEWMAT::ColumnVector grad, gradtmp, dist(2);
    NEWMAT::Matrix Disp, Dispinv, Data, weights;

    gradtmp.ReSize(2); gradtmp = 0.0;
    Disp.ReSize(querypoints.size(),2); Disp = 0;
    weights.ReSize(querypoints.size(),querypoints.size()); weights = 0;
    Data.ReSize(querypoints.size(),1);

    newresampler::Point point = SOURCE.get_coord(index);
    newresampler::Point origin = tangs.e1 * tangs.e2;
    origin.normalize(); origin *= RAD;

    newresampler::project_point(point - origin, tangs, y11, y21);

    for(unsigned int i = 0; i < querypoints.size(); i++)
    {
        newresampler::Point cr = TARGET.get_coord(querypoints[i]);
        newresampler::project_point(cr - origin, tangs, x11, x21);
        dist(1) = x11 - y11;
        dist(2) = x21 - y21;

        if((dist(1) * dist(1) + dist(2) * dist(2))>0)
        {
            Disp(i+1,1) = x11;
            Disp(i+1,2) = x21;
            //THE DISTANCES SHOULD BE CALCULATED FOR THE TANGENT PLANE
            // we ideally want an analytica expression
            weights(i+1,i+1) =
                    refweight(1, querypoints[i] + 1) * exp(-(dist(1) * dist(1) + dist(2) * dist(2)) / (2 * min_sigma * min_sigma));
            SUM += weights(i+1,i+1);
            JPsim += sim.Peek(querypoints[i]+1,index+1) * weights(i+1,i+1);
        }
    }

    if(SUM > 0) JPsim /= SUM;

    for (unsigned int i = 0; i < querypoints.size(); i++)
        Data(i + 1, 1) = sim.Peek(querypoints[i] + 1, index + 1) - JPsim;
        // difference of query point to mean similarity

    //estimates fits a gradient to the set samples and their similarity values using WLS
    NEWMAT::Matrix tmp = ((Disp.t() * weights) * Disp);
    if(tmp.Determinant() > 1e-5)
    {
        Dispinv = (tmp.i() * Disp.t()) * weights;
        gradtmp = Dispinv * Data;
    }

    grad.ReSize(3); grad = 0.0;
    grad(1) += (gradtmp(1)*(tangs.e1).X + gradtmp(2) * (tangs.e2).X);
    grad(2) += (gradtmp(1)*(tangs.e1).Y + gradtmp(2) * (tangs.e2).Y);
    grad(3) += (gradtmp(1)*(tangs.e1).Z + gradtmp(2) * (tangs.e2).Z);

    current_sim(index+1) = JPsim;

    return grad;
}

void Rigid_cost_function::update_similarity() {

    sim.set_input(FEAT->get_input_data());
    sim.set_reference(FEAT->get_reference_data());

#pragma omp parallel for
    for (int i = 0; i < SOURCE.nvertices(); i++)
        sim.calculate_sim_column_nbh(i);
}

NEWMAT::ColumnVector Rigid_cost_function::Evaluate_SIMGradient(int i, const newresampler::Tangs& tangs) {

    MISCMATHS::SpMat<int> found(TARGET.nvertices(),1);
    newresampler::Point point = SOURCE.get_coord(i);
    vector<int> querypoints;
    bool update = false;
    NEWMAT::ColumnVector vecnew;

    if(nbh->nrows(i))
    {
        int t_ind = (*nbh)(i, 0);
        vecnew.ReSize(3); vecnew = 0;

        do {
            querypoints.clear();
            newresampler::Triangle closest_triangle = targettree->get_closest_triangle(point);
            int n0 = closest_triangle.get_vertex_no(0),
                n1 = closest_triangle.get_vertex_no(1),
                n2 = closest_triangle.get_vertex_no(2);

            if((refweight(1, n0 + 1) > 0) && (refweight(1, n1 + 1) > 0) && (refweight(1, n2 + 1) > 0))
            {
                if(get_all_neighbours(i, querypoints, point, n0, TARGET, nbh, found)) update = true;
                if(get_all_neighbours(i, querypoints, point, n1, TARGET, nbh, found) || update) update = true;
                if(get_all_neighbours(i, querypoints, point, n2, TARGET, nbh, found) || update) update = true;
                if (update)
                {
                    nbh->at(i) = querypoints;
                    sim.calculate_sim_column_nbh(i);
                }
                // calculates similarity gradient
                vecnew = WLS_simgradient(tangs, i, querypoints);
                vecnew = inweight(1, i + 1) * refweight(1, t_ind + 1) * vecnew;
                update = false;
            }
            else
            {
                update = false;
                current_sim(i+1) = 0;
            }
        } while(update);
    }

    return vecnew;
}

void Rigid_cost_function::rotate_in_mesh(double a1, double a2, double a3){

#pragma omp parallel for
    for (int index = 0; index < SOURCE.nvertices(); index++)
    {
        newresampler::Point cii = SOURCE.get_coord(index);
        NEWMAT::ColumnVector V(3), VR(3);
        V(1) = cii.X; V(2) = cii.Y; V(3) = cii.Z;
        VR = newresampler::euler_rotate(V, a1, a2, a3);
        newresampler::Point ppc(VR(1), VR(2), VR(3));
        SOURCE.set_coord(index,ppc);
    }
}

double Rigid_cost_function::rigid_cost_mesh(double dw1, double dw2, double dw3){

    double SUM = 0;
    newresampler::Mesh tmp = SOURCE;

    rotate_in_mesh(dw1, dw2, dw3);

#pragma omp parallel for
    for (int index = 0; index < SOURCE.nvertices(); index++)
    {
        newresampler::Tangs T = calculate_tangs(index, SOURCE);
        Evaluate_SIMGradient(index,T);
#pragma omp critical
        SUM += current_sim(index + 1);
    }

    SOURCE = tmp;
    return SUM;
}

newresampler::Mesh Rigid_cost_function::run(){

    if(verbosity) std::cout << "Rigid registration started" << std::endl;

    double Euler1 = 0.0, Euler2 = 0.0, Euler3 = 0.0,
        mingrad_zero = 0.0, per = 0.0, RECinit = 0.0, RECfinal = 0.0;
    int min_iter = 0, loop = 0;

    NEWMAT::Matrix Rot(3,3), Rot_f(3,3);
    Rot_f = 0; Rot_f(1,1) = 1; Rot_f(2,2) = 1; Rot_f(3,3) = 1;
    NEWMAT::ColumnVector E(3);
    newresampler::Mesh tmp;

    double grad_zero = rigid_cost_mesh(Euler1, Euler2, Euler3);
    // tries different small rotations and gets mesh similarity for each of these points

    mingrad_zero = grad_zero;
    RECinit = grad_zero;

    while (spacing > 0.05)
    {
        double step = stepsize;
        per = spacing;

        for (int _iter = 1; _iter <= iters; _iter++)
        {
            //Initially Euler angles are zero
            Euler1 = 0.0; Euler2 = 0.0; Euler3 = 0.0;

            newresampler::Point grad;
            grad.X = (rigid_cost_mesh(Euler1 + per, Euler2, Euler3) - grad_zero) / per;
            grad.Y = (rigid_cost_mesh(Euler1, Euler2 + per, Euler3) - grad_zero) / per;
            grad.Z = (rigid_cost_mesh(Euler1, Euler2, Euler3 + per) - grad_zero) / per;
            grad.normalize();

            Euler1 = Euler1 + step * grad.X;
            Euler2 = Euler2 + step * grad.Y;
            Euler3 = Euler3 + step * grad.Z;

            if(verbosity)
            {
                cout << "******:   " << _iter << " per " << per << " loop" << loop
                     << " (loop*iters)+ _iter " << (loop * iters) + _iter << '\n'
                     << " grad_zero: " << grad_zero << '\n' << " mingrad_zero: " << mingrad_zero
                     << " min_iter " << min_iter << '\n' << " stepsize: " << step << endl;
            }

            tmp = SOURCE;

            rotate_in_mesh(Euler1, Euler2, Euler3);
            grad_zero = rigid_cost_mesh(Euler1, Euler2, Euler3);

            if(grad_zero > mingrad_zero)
            {
                mingrad_zero = grad_zero;
                min_iter = (loop * iters) + _iter;
                RECfinal = mingrad_zero;
            }

            if((loop * iters) + _iter - min_iter > 0)
            {
                step *= 0.5;
                SOURCE = tmp;
            }

            if(step < 1e-3) break;
        }
        loop++;
        spacing *= 0.5;
    }

    E(1) = Euler1;
    E(2) = Euler2;
    E(3) = Euler3;

    if(verbosity && (RECfinal > 0))
        std::cout << " Affine improvement: " << abs(((RECfinal-RECinit))/RECfinal*100) << "%" << std::endl;

    return SOURCE;
}

} //namespace newmeshreg
