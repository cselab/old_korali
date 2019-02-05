#include "coupled_ode_system.hpp"
#include "system_utils.hpp"

#include <boost/numeric/odeint.hpp>
#include <boost/bind.hpp>
#include <boost/array.hpp>


void CoupledOdeSystem::setObservations (const vec_d & times, 
    const std::vector<vec_d> & observations) {

    _ntimes = times.size();
    _times  = vec_d(_ntimes);
    _obsdim = observations.size();
    _obs    = std::vector<vec_d>(_obsdim);
    _sim    = std::vector<vec_d>(0);


    for(int i = 0; i < _ntimes; ++i) { _times[i] = times[i]; };

    for(int i = 0; i < _obsdim; ++i) {
        _obs[i]   = vec_d(_ntimes);
        for(int j = 0; j < _ntimes; ++j) _obs[i][j] = observations[i][j];
    }

}


vec_s CoupledOdeSystem::getIC() const
{
    vec_s ic = getModelIC();
    ic.resize(_dim + _dim * _numparam);
    try {
        vec_d grads;
        for(size_t i = 0; i < _dim; ++i) {
            stan::math::set_zero_all_adjoints();
            ic[i].grad();
            for(size_t j = 0; j < _numparam; ++j) 
                ic[_dim + _numparam * i + j] = _params[j].adj();
        }
    } catch (const std::exception& e) {
        std::cout << "coupled_ode_system.cpp::getIC() exception caught: " << e.what() << std::endl;
    }

    return ic;
}

void CoupledOdeSystem::step(const vec_d & z, vec_d & dz, double t) 
{

    /*
    vec_s z_s  = vec_s(_dim);
    vec_s dz_s = vec_s(_dim);
    for(int i = 0; i < _dim; ++i) {
        z_s[i]  = z[i];  //TODO:optimize? (DW)
        dz_s[i] = dz[i];
    } */
    vec_s z_s(z.begin(), z.end());
    vec_s dz_s(dz.begin(), dz.end());
    this->step(z_s, dz_s, t);
    
    for(int i = 0; i< dz_s.size(); ++i) dz[i] = value_of(dz_s[i]);


}

void CoupledOdeSystem::step(const vec_s & z, vec_s & dzOut, double t)
{
    const int COUPLED_DIM = _dim * (_numparam + 1);
    printf("t %lf\n",t);
    //printvec_s("z",z);

    evalModel(dzOut, z, t); //i think this is not used (DW)
    dzOut.resize(COUPLED_DIM, 0.0 );
    
    stan::math::start_nested();
    
    vec_d z_d  = vec_d(COUPLED_DIM);
    vec_d dz_d = vec_d(COUPLED_DIM) ;
    for(int i = 0; i < _dim; ++i) {
        z_d[i]  = value_of(z[i]);  //TODO:optimize? (DW)
        dz_d[i] = value_of(dzOut[i]);
    }
    
    Eigen::Map<const Eigen::Matrix<double,-1,-1, Eigen::ColMajor>>
            S_trans(&z_d[_dim], _numparam, _dim);
    Eigen::Map<Eigen::Matrix<double,-1,-1, Eigen::ColMajor>>
            GK_trans(&dz_d[_dim], _numparam, _dim);

    vec_s tmp = _params;
    for(int i = 0; i < _numparam; ++i) _params[i] = tmp[i].val();

    vec_s model_dot(_dim);

    vec_s z_in(z.begin(), z.end());

    //printvec_s("z_in",z_in);
    evalModel(model_dot, z_in, t);
    for(size_t i = 0; i < _dim; ++i) {
        stan::math::set_zero_all_adjoints_nested();

        model_dot[i].grad();

        for(size_t j = 0; j < _numparam; ++j) {
            _B_temp_trans(j,i) = _params[j].adj();
            //printf("_B_temp_trans(%zu,%zu) = %lf\n", j, i, _B_temp_trans(j,i));
        }

        for(size_t k = 0; k < _dim; ++k) {
            _A_trans(k,i) = z_in[k].adj();
            //printf("_A_trans(%zu,%zu) = %lf\n", k, i, _A_trans(k,i));
        }
    }

    stan::math::recover_memory_nested();
    
    GK_trans = S_trans*_A_trans+_B_temp_trans;
    
    for(int i = _dim; i< COUPLED_DIM; ++i) dzOut[i] = dz_d[i];
    printvec_s("dzOut (out)", dzOut);
    
}

return_type * CoupledOdeSystem::fitfun(double *x, int n, void* output, int *info)
{

    const int indexSigma = n-1;

    const double sigma     = x[indexSigma];
    const double sigma2    = sigma*sigma;
    const double invsigma3 = 1.0/(sigma2*sigma);

    vec_d theta_d(x, x+n);
    vec_s theta_s(x, x+n);

    setParams(theta_d);

    stan::math::start_nested();
    std::vector<vec_s> observable_s(_obsdim);

    vec_s ic_s = getIC();
    vec_d ic(ic_s.size());

    for(size_t i = 0; i < ic_s.size(); ++i) ic[i] = value_of(ic_s[i]);

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

    std::vector<vec_s> eq_sol(_obsdim);
    vec_d gradients;
    //associate sensitivities with respective values and find observable (using mala)
    for(size_t i = 0; i < _obsdim; ++i) {
        eq_sol[i].resize( _dim );
        for(size_t j = 0; j < _dim; j++) {
            gradients = vec_d(sensitivity_d[i].begin()+j*_numparam,
                              sensitivity_d[i].begin()+(j+1)*_numparam);
            eq_sol[i][j] = stan::math::precomputed_gradients(
                               equation_solution_d[i][j], theta_s, gradients);
        }
        vec_s temp_obs; // = findObservable(eq_sol[i]);
        observable_s.push_back(temp_obs);
    }

    scalar_t llk = 0;
    for(size_t i = 0; i < _ntimes; ++i) {
        for(size_t j = 0; j < _obsdim; ++j) {
            llk += (_obs[i][j]-observable_s[i][j]) *
                   (_obs[i][j]-observable_s[i][j]);
        }
    }

    //--------------------------------------------------------------------------------------

    llk /= sigma2;
    llk += _obsdim*(std::log(2 * M_PI) + std::log(sigma2));
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
        grad_loglike[i] = theta_s[i].adj();
    }

    for (size_t i = 0; i < _numparam; ++i) {
        gsl_vector_set(gradient, i, grad_loglike[i] );
        if (!std::isfinite(gsl_vector_get(gradient, i))) grad_err = 1;
    }

    double totObs = _ntimes * _obsdim;
    double tmp    = -totObs/sigma + sumVal*invsigma3;
    gsl_vector_set(gradient, indexSigma, tmp);
    if (!std::isfinite(gsl_vector_get(gradient, indexSigma))) grad_err = 1;

    return_type* result = reinterpret_cast<return_type*>(calloc(1, sizeof(return_type)) );

    // find Fischer Matrix
    Eigen::MatrixXd eigS(_numparam, _obsdim);
    Eigen::MatrixXd eigFIM = Eigen::MatrixXd::Zero(n,n);

    for(size_t i = 0; i < _obsdim; ++i) {
        stan::math::set_zero_all_adjoints_nested();
        observable_s[i][0].grad();
        for(size_t j = 0; j < _numparam; ++j) {
            eigS(j,i) = theta_s[j].adj();
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
    const vec_d& y_in,                  /* intial condition */
    const double integration_dt,
    double relative_tolerance,
    double absolute_tolerance,
    int max_num_steps)
{
    using boost::numeric::odeint::make_controlled;
    using boost::numeric::odeint::runge_kutta_dopri5;
    using boost::numeric::odeint::max_step_checker;
    
    auto systemLambda = [this] (const vec_d & y, vec_d & dy, double t) 
                            { return this->step(y, dy, t); };
    
    auto observerLambda = [this] (const vec_d & x, double t) 
                            { return this->observer(x, t); };
   
    vec_d y_test = vec_d(y_in.size()); //TODO: rmv later (DW)
    for(int i = 0; i < y_in.size(); ++i) y_test[i] = y_in[i];

    try {
        boost::numeric::odeint::integrate_times(
            make_controlled( absolute_tolerance, relative_tolerance, runge_kutta_dopri5<vec_d>() ),
            systemLambda,                           
            y_test,
            std::begin(_times),                     
            std::end(_times),     
            integration_dt,                
            observerLambda,
            max_step_checker(max_num_steps)
        );
    } catch (std::exception& e) {
        return std::make_pair(_sim, false);
    }
    return std::make_pair(_sim,true);
}


