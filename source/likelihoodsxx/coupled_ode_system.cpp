#include "coupled_ode_system.hpp"
#include "fitfun.hpp"
#include "system_utils.hpp"

#include <boost/numeric/odeint.hpp>
#include <boost/bind.hpp>
#include <boost/array.hpp>

namespace fitfun {

void CoupledOdeSystem::setObservations (const vec_d & times,
                                        const std::vector<vec_d> & observations)
{

    _ntimes = times.size();
    _times  = vec_d(times.begin(), times.end());
    _obsdim = observations.begin()->size();
    _obs    = std::vector<vec_d>(_ntimes);

    for(size_t i = 0; i < _ntimes; ++i) {
        _obs[i]   = vec_d(_obsdim);
        for(size_t j = 0; j < _obsdim; ++j) _obs[i][j] = observations[i][j];
    }

}


vec_d CoupledOdeSystem::getIC(const vec_d & params) const
{
    vec_s params_s(params.begin(), params.end());
    vec_s ic_s = getModelIC_s(params_s);

    if(_mala) {
        ic_s.resize(_dim + _dim * _numparam);
        try {
            vec_d grads;
            for(size_t i = 0; i < _dim; ++i) {
                stan::math::set_zero_all_adjoints();
                ic_s[i].grad();
                for(size_t j = 0; j < _numparam; ++j)
                    ic_s[_dim + _numparam * i + j] = params_s[j].adj();
            }
        } catch (const std::exception& e) {
            printf("CoupledOdeSystem::getIC() exception caught: %s",e.what());
        }
    }

    vec_d ic(ic_s.size());
    for(int i = 0; i < ic_s.size(); ++i) ic[i] = ic_s[i].val();
    return ic;
}


void CoupledOdeSystem::observer(const vec_d & state, double t)
{
    vec_d sol(state.begin(), state.end());
    _sim.push_back(sol);
}


void CoupledOdeSystem::step(const vec_d & z, vec_d & dzOut, double t)
{
    //printf("t %lf\n",t);
    //printvec_d("z", z);

    if(_mala) stan::math::start_nested();

    vec_s params_s( _params.begin(), _params.end() );
    vec_s z_s( z.begin(), z.begin()+_dim );
    vec_s model_dot_s(_dim);

    evalModel_s(model_dot_s, z_s, params_s, t);
    for(size_t i =0 ; i< _dim; ++i) dzOut[i] = model_dot_s[i].val();

    if (_mala) {

        const int COUPLED_DIM = _dim * (_numparam + 1);
        dzOut.resize(COUPLED_DIM, 0.0 );


        Eigen::Map<const Eigen::Matrix<double,-1,-1, Eigen::ColMajor>>
                S_trans(&z[_dim], _numparam, _dim);
        Eigen::Map<Eigen::Matrix<double,-1,-1, Eigen::ColMajor>>
                GK_trans(&dzOut[_dim], _numparam, _dim);

        //printvec_s("z_s", z_s);

        for(size_t i = 0; i < _dim; ++i) {
            stan::math::set_zero_all_adjoints_nested();

            model_dot_s[i].grad();

            for(size_t j = 0; j < _numparam; ++j) {
                _B_temp_trans(j,i) = params_s[j].adj();
                //printf("_B_temp_trans(%zu,%zu) = %lf\n", j, i, _B_temp_trans(j,i));
            }

            for(size_t k = 0; k < _dim; ++k) {
                _A_trans(k,i) = z_s[k].adj();
                //printf("_A_trans(%zu,%zu) = %lf\n", k, i, _A_trans(k,i));
            }
        }

        stan::math::recover_memory_nested();

        GK_trans = S_trans*_A_trans+_B_temp_trans;
    }

    //printvec_d("dzOut", dzOut);
    return;
}

std::pair<std::vector<vec_d >, bool> CoupledOdeSystem::integrate_boost(
    vec_d & y_in,
    const double integration_dt,
    double relative_tolerance,
    double absolute_tolerance,
    int max_num_steps)
{
    using boost::numeric::odeint::make_controlled;
    using boost::numeric::odeint::runge_kutta_dopri5;
    using boost::numeric::odeint::max_step_checker;

    //TODO: can we work with plain lambdas without the step and observer
    //function? (DW) (and maybe get rid of _sim)
    auto systemLambda = [this] (const vec_d & y, vec_d & dy, double t) {
        return this->step(y, dy, t);
    };

    auto observerLambda = [this] (const vec_d & x, double t) {
        return this->observer(x, t);
    };

    int idx;
    vec_d t(2);
    _sim = std::vector<vec_d>(0);
    std::vector<vec_d> sol(_times.size());

    try {

        if(_t0 < _times[0]) {
            idx  = 0;
            t[0] = _t0;
            t[1] = _times[0];
        } else { 
            idx    = 1;
            t[0]   = _times[0];
            sol[0] = y_in;
        }
        
        for( ; idx < _times.size(); ++idx) {
    
            t[1] = _times[idx];
                   
            boost::numeric::odeint::integrate_times(
                make_controlled( absolute_tolerance, relative_tolerance, runge_kutta_dopri5<vec_d>() ),
                systemLambda,
                y_in,
                std::begin(t),
                std::end(t),
                integration_dt,
                observerLambda,
                max_step_checker(max_num_steps)
            );

            t[0]     = _times[idx];
            y_in     = _sim.back();
            sol[idx] = _sim.back();
        }

    } catch (std::exception& e) {
        printf("CoupledOdeSystem::integrate_boost: \
                Exception caught in integrate_times: %s", e.what());
        return std::make_pair(sol, false);
    }

    return std::make_pair(sol,true);
}


double CoupledOdeSystem::evaluate(const double *x, int n, void* output, int *info)
{
    const int indexSigma = n-1;

    const double sigma     = x[indexSigma];
    const double sigma2    = sigma*sigma;
    const double invsigma3 = 1.0/(sigma2*sigma);

    return_type* result = static_cast<return_type*>(output);

    vec_d params(x, x+indexSigma);
    setParams(params);

    stan::math::start_nested();

    vec_d ic = getIC(params);

    bool success;
    std::vector<vec_d > solution_of_coupled_system;
    std::tie(solution_of_coupled_system, success) = integrate_boost(ic);

    if(!success) {
        printf("CoupledOdeSystem::fitfun : error in integrate_boost(ic)\n");
        result->error_flg = 1;
        result->loglike   = -1e6;
        stan::math::recover_memory_nested();
        return result->loglike;
    }

    vec_s params_s(params.begin(), params.end());
    std::vector<vec_s> observable_s(_ntimes);
    std::vector<vec_d> equation_solution;
    std::vector<vec_d> equation_sensitivity;
    decouple(solution_of_coupled_system, equation_solution,
             equation_sensitivity, _dim, _numparam);


    std::vector<vec_s> eq_sol(_ntimes);
    vec_d gradients;
    for(size_t i = 0; i < _ntimes; ++i) {
        //printf("i: %zu\n", i); printvec_d("equation_solution", equation_solution[i]);
        if(_mala) {
            eq_sol[i]    = vec_s(_dim, 0.0);
            for(size_t j = 0; j < _dim; ++j) {
                gradients    = vec_d(equation_sensitivity[i].begin()+j*_numparam,
                                     equation_sensitivity[i].begin()+(j+1)*_numparam);
                eq_sol[i][j] = stan::math::precomputed_gradients(
                                   equation_solution[i][j], params_s, gradients);
            }
        } else {
            eq_sol[i] = vec_s(equation_solution[i].begin(), equation_solution[i].end());
        }
        observable_s[i] = calculateObservable(eq_sol[i]);
    }

    scalar_t llk = 0;
    for(size_t i = 0; i < _ntimes; ++i) {
        for(size_t j = 0; j < _obsdim; ++j) {
            llk += (_obs[i][j]-observable_s[i][j]) *
                   (_obs[i][j]-observable_s[i][j]);
        }
    }

    llk /= sigma2;
    llk += _obsdim*(std::log(2 * M_PI) + std::log(sigma2));
    llk *= -0.5;

    double llkVal = llk.val();

    if ( !std::isfinite(llkVal) ) {
        printf("CoupledOdeSystem::fitfun : log likelihood not finite\n");
        result->error_flg   = 1;
        result->loglike     = -1e5;
        return result->loglike;
    }

    result->loglike = llkVal;
    if(!_mala) {
        stan::math::recover_memory_nested();
        return result->loglike;
    }

    /* FOR MALA BELOW */

    int grad_err = 0;
    gsl_vector *gradient = gsl_vector_alloc(n);
    vec_d grad_loglike(_numparam);

    stan::math::set_zero_all_adjoints_nested();
    llk.grad();
    for(size_t i = 0; i < _numparam; ++i) {
        grad_loglike[i] = params_s[i].adj();
    }

    for (size_t i = 0; i < _numparam; ++i) {
        gsl_vector_set(gradient, i, grad_loglike[i] );
        if (!std::isfinite(gsl_vector_get(gradient, i))) grad_err = 1;
    }

    double totObs = _ntimes * _obsdim;
    double tmp    = -totObs/sigma + llkVal*invsigma3;
    gsl_vector_set(gradient, indexSigma, tmp);
    if (!std::isfinite(gsl_vector_get(gradient, indexSigma))) grad_err = 1;

    if (grad_err != 0) {
        printf("CoupledOdeSystem::fitfun: Gradient not finite. \n");
        result->error_flg = 1;
        result->loglike   = llkVal;
        gsl_vector_free(gradient);
        stan::math::recover_memory_nested();
        return result->loglike;
    }

    // find Fischer Matrix
    Eigen::MatrixXd eigS(_numparam, _obsdim);
    Eigen::MatrixXd eigFIM = Eigen::MatrixXd::Zero(n,n);

    for(size_t i = 0; i < _obsdim; ++i) {
        stan::math::set_zero_all_adjoints_nested();
        observable_s[i][0].grad();
        for(size_t j = 0; j < _numparam; ++j) {
            eigS(j,i) = grad_loglike[j];
        }
    }

    // S*S^T = S2
    eigFIM.block(0, 0, _numparam, _numparam) = eigS*eigS.transpose();
    eigFIM(indexSigma,indexSigma)  = 2*totObs;

    Eigen::MatrixXd I = Eigen::MatrixXd::Identity(n,n);
    Eigen::MatrixXd eigInv_FIM =
        eigFIM.fullPivHouseholderQr().solve(I)*sigma2; // TODO: check, we scale here and some lines below again?? (DW)

    // TODO: pass display through info array? (DW)
    // if (display>2) {
    // bool a_solution_exists = (eigFIM*eigInv_FIM).isApprox(I*sigma2, 1e-2);
    // double relative_error  = (eigFIM*eigInv_FIM - I*sigma2).norm() / I.norm(); // norm() is L2 norm
    // printf("CoupledOdeSystem::fitfun: eigFIM inversion succesfull: %d (rel error: %lf)\n", a_solution_exists, relative_error);
    // }
    
    gsl_matrix * inv_FIM  = gsl_matrix_calloc(n, n);
    for(int i = 0; i< n; ++i) {
        for(int j = 0; j<n; ++j) {
            gsl_matrix_set(inv_FIM, i, j, eigInv_FIM(i,j));
        }
    }

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigEs(n);
    eigEs.compute(eigFIM);
    Eigen::VectorXd eigInvEigenValues  = eigEs.eigenvalues().cwiseInverse()*sigma2; //scaling of matrix is scaling of eigenvalues
    Eigen::MatrixXd eigInvEigenVectors = eigEs.eigenvectors(); //A^(-1)^(-1) = A from FIM = ADA^(-1)

    gsl_matrix *evec = gsl_matrix_calloc(n, n);
    gsl_vector *eval = gsl_vector_calloc(n);

    for(int i = 0; i < n; ++i) {
        gsl_vector_set(eval, i, eigInvEigenValues(i));
        for(int j = 0; j < n; ++j) {
            gsl_matrix_set(evec, i, j, eigInvEigenVectors(i,j));
        }
    }

    bool posdef = (gsl_vector_min(eval) > 0.0);

    result->grad      = gradient;
    result->posdef    = posdef;
    result->error_flg = 0; // error_flg = false

    gsl_vector_free(gradient);

    for (int i = 0; i< n; ++i) {
        if (!std::isfinite(gsl_vector_get(eval, i))) {
            printf("CoupledOdeSystem::fitfun: Eval not finite. \n");
            result->error_flg = 2;
            stan::math::recover_memory_nested();
            gsl_matrix_free(inv_FIM);
            gsl_matrix_free(evec);
            gsl_vector_free(eval);
            return result->loglike;
        }
    }

    result->eval = eval;

    for (int i = 0; i<n; ++i) {
        for (int j = 0; j<n; ++j) {
            if (!std::isfinite(gsl_matrix_get(evec, i, j)) ||
                    !std::isfinite(gsl_matrix_get(inv_FIM, i, j)) )  {
                printf("CoupledOdeSystem::fitfun: Evec or inv_FIM not finite. \n");
                result->error_flg = 2;
                stan::math::recover_memory_nested();
                gsl_matrix_free(inv_FIM);
                gsl_matrix_free(evec);
                gsl_vector_free(eval);
                return result->loglike;
            }
        }
    }

    result->cov  = inv_FIM;
    result->evec = evec;

    stan::math::recover_memory_nested();
    gsl_matrix_free(inv_FIM);
    gsl_matrix_free(evec);
    gsl_vector_free(eval);

    return result->loglike;
}

}//namespace fitfun
