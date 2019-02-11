#ifndef _COUPLED_ODE_SYSTEM_HPP_
#define _COUPLED_ODE_SYSTEM_HPP_

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

#include "system_utils.hpp"
#include "Ifitfun.hpp"

namespace fitfun {

class CoupledOdeSystem : public IFitfun
{

public:

    CoupledOdeSystem(int numparam, int odedim, bool mala)
        : _obsdim(0), _numparam(numparam),  _dim(odedim),  _mala(mala)
    {
        _A_trans      = Eigen::MatrixXd::Zero(_dim, _dim);
        _B_temp_trans = Eigen::MatrixXd::Zero(_numparam, _dim);
    }

    double evaluate(const double *x, int n, void* output, int *info);
    void initialize(int argc, const  char **argv) {};
    void finalize() {};

    void setParams(vec_d params)
    {
        _params = params;
    };
    void setObservations (const vec_d & times, const std::vector<vec_d> & observations);
    void step(const vec_d & z, vec_d & dz, double t);

    vec_d getIC(const vec_d & params) const;

protected:

    int _obsdim;

    std::vector<vec_d> _obs; // [_obsdim] x [_ntimes]
    std::vector<vec_d> _sim; // [_dim] x [_ntimes]

    virtual vec_s getModelIC_s(const vec_s & params) const = 0;
    virtual vec_s calculateObservable(const vec_s & solution) const = 0;
    virtual void evalModel_s(vec_s & dyOut, const vec_s & y, const vec_s & params, double t) = 0;


private:

    const int _numparam;
    const int _dim;
    
    int _ntimes;
    vec_d _times;

    vec_d _params;
    const bool _mala;

    Eigen::MatrixXd _A_trans;
    Eigen::MatrixXd _B_temp_trans;

    void observer(const vec_d & state, double t);

    std::pair<std::vector<vec_d >, bool> integrate_boost(
        vec_d & y_in,
        const double integration_dt = 0.1,
        double relative_tolerance = 1e-8,
        double absolute_tolerance = 1e-8,
        int max_num_steps = 1e4);

};

}//namespace fitfun

#endif// _COUPLED_ODE_SYSTEM_HPP_