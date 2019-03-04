// libgp - Gaussian process library for Machine Learning
// Copyright (c) 2013, Manuel Blum <mblum@informatik.uni-freiburg.de>
// All rights reserved.

#include "cov_periodic.h"
#include <cmath>

namespace libgp
{

  CovPeriodic::CovPeriodic() {}

  CovPeriodic::~CovPeriodic() {}

  bool CovPeriodic::init(int n)
  {
    input_dim = n;
    param_dim = 3;
    loghyper.resize(param_dim);
    loghyper.setZero();
    return true;
  }

  double CovPeriodic::get(const Eigen::VectorXd &x1, const Eigen::VectorXd &x2)
  {
    double s = sin(M_PI * (x1-x2).norm() / T) / ell;
    return sf2*exp(-2*s*s);
  }

  void CovPeriodic::grad(const Eigen::VectorXd &x1, const Eigen::VectorXd &x2, Eigen::VectorXd &grad)
  {
    double k = M_PI * (x1-x2).norm() / T;
    double s = sin(k) / ell;
    // why is the last elements 0?
    // the commented derivative corresponds to derivative w.r.t. log(T)
    grad << 4*sf2*exp(-2*s*s)*s*s, 2*sf2*exp(-2*s*s), 0;// 4*sf2/ell*exp(-2*s*s)*s*cos(k)*k;
    // grad << 4*sf2*exp(-2*s*s)*s*s, 2*sf2*exp(-2*s*s), 4*sf2/ell*exp(-2*s*s)*s*cos(k)*k;

  }


  void CovPeriodic::gradx(const Eigen::VectorXd &x1, const Eigen::VectorXd &x2, Eigen::VectorXd &grad)
  {
    double c = 2*ell*M_PI/T;
    double z = M_PI * (x1-x2).norm() / T ;
    double f = sf2*exp(-2*z*z);
    double k = c*c*f*sin(z)*cos(z)/z;
    grad << k*(x2-x1);
  }



  void CovPeriodic::set_loghyper(const Eigen::VectorXd &p)
  {
    CovarianceFunction::set_loghyper(p);
    ell = exp(loghyper(0));
    sf2 = exp(2*loghyper(1));
    T = fabs(loghyper(2));
  }

  std::string CovPeriodic::to_string()
  {
    return "CovPeriodic";
  }

}
