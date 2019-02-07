#ifndef _COUPLED_ODE_SYSTEM_HPP_
#define _COUPLED_ODE_SYSTEM_HPP_

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

#include <stan/math.hpp>

#include "fitfun.hpp"

typedef stan::math::var scalar_t;

typedef std::vector<double>  vec_d;
typedef std::vector<scalar_t>  vec_s;

class CoupledOdeSystem : public Fitfun
{

public:

    CoupledOdeSystem(int numparam, int odedim, bool mala)
        : _numparam(numparam),  _dim(odedim), _obsdim(0), _it(0), _mala(mala)
    {
        _A_trans      = Eigen::MatrixXd::Zero(_dim, _dim);
        _B_temp_trans = Eigen::MatrixXd::Zero(_numparam, _dim);
    }

    return_type * fitfun(double *x, int n, void* output, int *info);
   
    void setParams(vec_d params) { _params = params; };
    void setObservations (const vec_d & times, const std::vector<vec_d> & observations);
    void setStartTime (double it) { _it = it; };
    void step(const vec_s & z, vec_s & dz, double t);
    void step(const vec_d & z, vec_d & dz, double t);


    vec_d getIC(const vec_d & params) const;
    vec_s getIC(const vec_s & params) const;

protected:

    const int _numparam;
    const int _dim;
    int _obsdim;
    double _it;
 
    std::vector<vec_d> _obs; // [_obsdim] x [_ntimes]
    std::vector<vec_d> _sim; // [_dim] x [_ntimes]
    
    inline int getDim() const { return _dim; };

    virtual vec_s getModelIC_s(const vec_s & params) const = 0;
    virtual void evalModel_s(vec_s & dyOut, const vec_s & y, const vec_s & params, double t) = 0;

    virtual vec_d getModelIC(const vec_d & params) const = 0;
    virtual void evalModel(vec_d & dyOut, const vec_d & y, const vec_d & params, double t) = 0;

    virtual vec_s calculateObservable(const vec_s & solution) = 0;

private:

    int _ntimes;    
    vec_d _times;
    
    vec_d _params;
    bool _mala;
    
    Eigen::MatrixXd _A_trans;
    Eigen::MatrixXd _B_temp_trans;

    void observer(const vec_d & state, double t);
    
    std::pair<std::vector<vec_d >, bool> integrate_boost(
        const vec_d& y_in,
        const double integration_dt = 0.1,
        double relative_tolerance = 1e-8,
        double absolute_tolerance = 1e-8,
        int max_num_steps = 1e4);

};


#endif// _COUPLED_ODE_SYSTEM_HPP_
