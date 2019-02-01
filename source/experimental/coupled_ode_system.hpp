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

    CoupledOdeSystem(int numparam, int odedim, double it)
        : _numparam(numparam), _obsdim(0), _dim(odedim), _it(it)
    {
        _trans     = Eigen::MatrixXd::Zero(_dim, numparam);
        _trans_tmp = Eigen::MatrixXd::Zero(_dim, numparam);
    }

    return_type * fitfun(double *x, int n, void* output, int *info);
    
    void setObservations (const vec_d & times, const std::vector<vec_d> & observations);
    void operator() (const vec_s & z, vec_s & dz, double t);
    void operator() (const vec_d & z, vec_d & dz, double t);

protected:

    std::vector<vec_d> _obs; // [_obdsdim] x [_ntimes]
    std::vector<vec_d> _sim; // [_obdsdim] x [_ntimes]

    virtual vec_s getModelIC(const vec_s & params) const = 0;
    virtual void evalModel(vec_s & dyOut, const vec_s & y, double t) = 0;

    inline int getDim() const { return _dim; };

private:

    int _numparam;
    int _ntimes;    
    int _obsdim;
    int _dim;
    double _it;
    
    vec_d _times;
    
    Eigen::MatrixXd _trans;
    Eigen::MatrixXd _trans_tmp;

    vec_s getIC(const vec_s & params) const;
    
    std::pair<std::vector<vec_d >, bool> integrate_boost(
        const vec_d& y_in,
        const double integration_dt = 0.1,
        double relative_tolerance = 1e-8,
        double absolute_tolerance = 1e-8,
        int max_num_steps = 1e3);

    void observer(const vec_d & coupled_state, double t);
};


#endif// _COUPLED_ODE_SYSTEM_HPP_
