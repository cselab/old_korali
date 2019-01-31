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

    CoupledOdeSystem(int dim, int numparam, vec_d params, double it, int nobs)
        : _dim(dim), _numparam(numparam), _it(it), _numobs(nobs)
    {
        _params = vec_s(_numparam);
        for(int i = 0; i< _numparam; ++i) {
            _params[i] = params[i];
        };
        _trans     = Eigen::MatrixXd::Zero(_dim, _numparam);
        _trans_tmp = Eigen::MatrixXd::Zero(_dim, _numparam);
    }

    vec_s getIC(const vec_s & params) const;

    return_type * fitfun(double *x, int n, void* output, int *info);

    void operator() (const vec_s & z, vec_s & dz, double t);

protected:

    int _dim;
    int _numparam;
    vec_s _params;

    double _it;

    vec_d _times;

    int _numobs;// dim of outer vec (num observations)
    std::vector<vec_d> _obs;

    virtual vec_s getModelIC(const vec_s & params) const = 0;
    virtual void evalModel(vec_s & dyOut, double t, const vec_s & y) = 0;

    inline int getDim() const
    {
        return _dim;
    };

private:

    Eigen::MatrixXd _trans;
    Eigen::MatrixXd _trans_tmp;

    std::pair<std::vector<vec_d >, bool> integrate_boost(
        const vec_d& y_in,
        double relative_tolerance = 1e-8,
        double absolute_tolerance = 1e-8,
        int max_num_steps = 1e3);

};


#endif// _COUPLED_ODE_SYSTEM_HPP_
