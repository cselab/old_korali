#include "coupled_ode_system.hpp"
#include "system_utils.hpp"

#include <boost/numeric/odeint.hpp>

vec_s CoupledOdeSystem::getIC(const vec_s & params) const
{

    vec_s ic = getModelIC(params);

    try {
        vec_d grads;
        for(size_t i = 0; i < _dim; ++i) {
            stan::math::set_zero_all_adjoints();
            ic[i].grad();
            for(size_t j = 0; j < _numparam; ++j) {
                ic.push_back( params[j].adj() );
            }
        }
    } catch (const std::exception& e) {
        std::cout << "coupled_ode_system.cpp: exception caught: " << e.what() << std::endl;
    }

    return ic;
}

void CoupledOdeSystem::operator() (const vec_s & z, vec_s & dz, double t)
{
    evalModel(dz, t, z);

    vec_d z_d  = vec_d(_dim);
    vec_d dz_d = vec_d(_dim);
    for(int i = 0; i < _dim; ++i) {
        z_d[i]  = value_of(z[i]);  //TODO:optimize? (DW)
        dz_d[i] = value_of(dz[i]);
    }
    vec_s z_in( z.begin(), z.begin()+_dim );
    vec_s model_dot( _dim );

    evalModel(model_dot, t, z);


    Eigen::Map<const Eigen::Matrix<double,-1,-1, Eigen::ColMajor>>
            S_trans(&z_d[_dim], _numparam, _dim);
    Eigen::Map<Eigen::Matrix<double,-1,-1, Eigen::ColMajor>>
            B_trans(&dz_d[_dim], _numparam, _dim);

    for(size_t i = 0; i < _dim; ++i) {
        stan::math::set_zero_all_adjoints_nested();

        model_dot[i].grad();
        for(size_t j = 0; j < _dim; ++j) {
            _trans_tmp(j,i) = _params[j].adj();
        }

        for(size_t k = 0; k < _dim; ++k) {
            _trans(k,i) = z[k].adj();
        }
    }

    stan::math::recover_memory_nested();

    B_trans = S_trans*_trans+_trans_tmp; //TODO: why calculate this? (DW)

}

return_type * CoupledOdeSystem::fitfun(double *x, int n, void* output, int *info)
{

    int indexSigma = n-1;

    double sigma     = x[indexSigma];
    double sigma2    = sigma*sigma;
    double invsigma3 = 1.0/(sigma2*sigma);

    std::vector<scalar_t> theta_s;
    for (size_t i = 0; i < n; i++) theta_s.push_back(x[i]);

    stan::math::start_nested();
    std::vector<vec_s> observable_s;
    observable_s.reserve(_numobs);


    vec_s ic_tmp = getIC(theta_s);
    vec_d ic;

    for(size_t i = 0; i < _dim; ++i) {
        ic.push_back(value_of(ic_tmp[i]));
    }

    bool success; //tmp (DW)
    std::vector<vec_d > solution_of_coupled_system;
    std::tie(solution_of_coupled_system, success) = integrate_boost(ic);

    if(!success) {
        return_type* result = reinterpret_cast<return_type*>(calloc(1, sizeof(return_type)) );
        result->error_flg = 1;
        result->loglike   = -1e6;
        stan::math::recover_memory_nested();
        return result;
    }

    std::vector<vec_d> equation_solution_d;
    std::vector<vec_d> sensitivity_d;
    decouple(solution_of_coupled_system, sensitivity_d, equation_solution_d,
             _dim, _numparam);

#ifdef OBSERVE_SENS
    vec_fitfun1 = vec_d(solution_of_coupled_system.back().begin(), solution_of_coupled_system.back().end());
#endif

    std::vector<vec_s> eq_sol(_numobs);
    vec_d gradients;
    //associate sensitivities with respective values and find observable (using mala)
    for(size_t i = 0; i < _numobs; ++i) {
        eq_sol[i].resize( _dim );
        for(size_t j = 0; j < _dim; j++) {
            gradients = vec_d(sensitivity_d[i].begin()+j*_numparam,
                              sensitivity_d[i].begin()+(j+1)*_numparam);
            eq_sol[i][j] = stan::math::precomputed_gradients(
                               equation_solution_d[i][j], _params, gradients);
        }
        vec_s temp_obs; // = findObservable(eq_sol[i]);
        observable_s.push_back(temp_obs);
    }

    scalar_t llk = 0;
    for(size_t i = 0; i < _times.size(); ++i) {
        for(size_t j = 0; j < _numobs; ++j) {
            llk += (_obs[i][j]-observable_s[i][j]) *
                   (_obs[i][j]-observable_s[i][j]);
        }
    }

    //--------------------------------------------------------------------------------------

    llk /= sigma2;
    llk += _numobs*(std::log(2 * M_PI) + std::log(sigma2));
    llk *= -0.5;

    double sumVal = llk.val();

    if ( !std::isfinite(sumVal) ) {
        return_type* result = reinterpret_cast<return_type*>(calloc(1, sizeof(return_type)) );
        result->error_flg   = 1;
        result->loglike     = -1e5;
        stan::math::recover_memory_nested();
        return result;
    }

    int grad_err = 0;
    gsl_vector *gradient = gsl_vector_alloc(n);
    vec_d grad_loglike(_numparam);

    stan::math::set_zero_all_adjoints_nested();
    llk.grad();
    for(size_t i = 0; i < _numparam; ++i) {
        grad_loglike[i] = 0;
        grad_loglike[i] = _params[i].adj();
    }

    for (size_t i = 0; i < _numparam; ++i) {
        gsl_vector_set(gradient, i, grad_loglike[i] );
        if (!std::isfinite(gsl_vector_get(gradient, i))) grad_err = 1;
    }

    double totObs = _times.size() * _numobs;
    double tmp    = -totObs/sigma + sumVal*invsigma3;
    gsl_vector_set(gradient, indexSigma, tmp);
    if (!std::isfinite(gsl_vector_get(gradient, indexSigma))) grad_err = 1;

    return_type* result = reinterpret_cast<return_type*>(calloc(1, sizeof(return_type)) );

    // find Fischer Matrix
    Eigen::MatrixXd eigS(_numparam, _numobs);
    Eigen::MatrixXd eigFIM = Eigen::MatrixXd::Zero(n,n);

    for(size_t i = 0; i < _numobs; ++i) {
        stan::math::set_zero_all_adjoints_nested();
        observable_s[i][0].grad();
        for(size_t j = 0; j < _numparam; ++j) {
            eigS(j,i) = _params[j].adj();
        }
    }

    // S*S^T = S2
    eigFIM.block(0, 0, _numparam, _numparam) = eigS*eigS.transpose();
    eigFIM(indexSigma,indexSigma)  = 2*totObs;

    if (grad_err != 0) {
        result->error_flg = 1;
        result->loglike   = sumVal;
        printf("Fitfun says: Error in FIM or grad. \n");
        gsl_vector_free(gradient);
        stan::math::recover_memory_nested();
        return result;
    }

    Eigen::MatrixXd eigInv_FIM =
        eigFIM.fullPivHouseholderQr().solve(Eigen::MatrixXd::Identity(n,n))*sigma2;

    gsl_matrix * inv_FIM  = gsl_matrix_calloc(n, n);
    for(std::size_t i = 0; i< n; ++i) {
        for(std::size_t j = 0; j<n; ++j) {
            gsl_matrix_set(inv_FIM, i, j, eigInv_FIM(i,j));
        }
    }

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigEs(n);
    eigEs.compute(eigFIM);
    Eigen::VectorXd eigInvEigenValues  = eigEs.eigenvalues().cwiseInverse()*sigma2; //scaling of matrix is scaling of eigenvalues
    Eigen::MatrixXd eigInvEigenVectors = eigEs.eigenvectors(); //A^(-1)^(-1) = A from FIM = ADA^(-1)

    gsl_matrix *evec = gsl_matrix_calloc(n, n);
    gsl_vector *eval = gsl_vector_calloc(n);

    for(std::size_t i = 0; i < n; ++i) {
        gsl_vector_set(eval,i, eigInvEigenValues(i));
        for(std::size_t j = 0; j < n; ++j) {
            gsl_matrix_set(evec, i, j, eigInvEigenVectors(i,j));
        }
    }

    bool posdef = (gsl_vector_min(eval) > 0.0);

    result->loglike   = sumVal;
    result->grad      = gradient;
    result->posdef    = posdef;
    result->error_flg = 0; // error_flg = false

    gsl_vector_free(gradient);

    for (int i = 0; i< n; ++i) {
        if (!std::isfinite(gsl_vector_get(eval, i))) {
            result->error_flg = 2;
            stan::math::recover_memory_nested();
            gsl_vector_free(eval);
            return result;
        }
    }

    result->eval = eval;

    for (int i = 0; i<n; ++i) {
        for (int j = 0; j<n; ++j) {
            if (!std::isfinite(gsl_matrix_get(evec, i, j)) ||
                    !std::isfinite(gsl_matrix_get(inv_FIM, i, j)) )  {
                result->error_flg = 2;
                stan::math::recover_memory_nested();
                gsl_matrix_free(inv_FIM);
                gsl_matrix_free(evec);
                gsl_vector_free(eval);
                return result;
            }
        }
    }

    result->cov  = inv_FIM;
    result->evec = evec;

    stan::math::recover_memory_nested();
    gsl_matrix_free(inv_FIM);
    gsl_matrix_free(evec);
    gsl_vector_free(eval);

    return result;
}


std::pair<std::vector<vec_d >, bool> CoupledOdeSystem::integrate_boost(
    const vec_d& y_in,
    double relative_tolerance,
    double absolute_tolerance,
    int max_num_steps)
{

    using boost::numeric::odeint::integrate_adaptive;
    using boost::numeric::odeint::integrate_times;
    using boost::numeric::odeint::make_controlled;
    using boost::numeric::odeint::runge_kutta_dopri5;
    using boost::numeric::odeint::max_step_checker;

    std::vector<vec_d > y_coupled(_times.size());
    Observer observer(y_coupled);

    const double step_size = 0.1;
    vec_d numerical_sol;
    for(size_t i = 0; i< y_in.size(); i++) {
        numerical_sol.push_back(y_in[i]);
    }

    observer(numerical_sol, _times[0]);

    int restart_idx = 0;
    double start_time, end_time;
    bool restart = false;
    bool callobserver = true;

    end_time = _times[0];
    vec_d integration_times(2);

    for(size_t i = 1; i < _times?gc.size(); i++) {
        start_time = end_time;
        bool no_skip = restart_idx < restart_points.size();
        if (no_skip && std::abs(restart_points[restart_idx]-_times?gc[i]) < 1e-8) {
            end_time = restart_points[restart_idx];
            restart = true;
            restart_idx++;
        } else if (no_skip && restart_points[restart_idx] < _times?gc[i] ) {
            end_time = restart_points[restart_idx];
            restart = true;
            callobserver = false;
            restart_idx++;
            if(i>0) i--;
        } else {
            end_time = _times?gc[i];
        }

        std::vector<vec_d> pseudo_vec(2*_times?gc.size()+2*restart_points.size());
        Observer pseudo_obs(pseudo_vec);

        integration_times[0] = start_time;
        integration_times[1] = end_time;

        try {
            // integrate_adaptive(
            //     make_controlled( absolute_tolerance, relative_tolerance, runge_kutta_dopri5<vec_d>() ),
            //     system,
            //     numerical_sol,
            //     start_time,
            //     end_time,
            //     step_size );
            integrate_times(
                make_controlled( absolute_tolerance, relative_tolerance, runge_kutta_dopri5<vec_d>() ),
                system,
                numerical_sol,
                boost::begin(integration_times),
                boost::end(integration_times),
                step_size,
                pseudo_obs,
                max_step_checker(max_num_steps)
            );
        } catch (std::exception& e) {
            return std::make_pair(y_coupled, false);
        }


        if(restart) {
            vec_d cpy_numerical_sol(numerical_sol.begin(), numerical_sol.end() );
            system.restart( cpy_numerical_sol, numerical_sol, end_time);
            restart = false;
        }
        if (callobserver) {
            observer(numerical_sol, end_time);
        } else {
            callobserver = true;
        }
    }

    return std::make_pair(y_coupled,true);
}


