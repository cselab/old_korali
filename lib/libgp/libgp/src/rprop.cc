// libgp - Gaussian process library for Machine Learning
// Copyright (c) 2013, Manuel Blum <mblum@informatik.uni-freiburg.de>
// All rights reserved.

#include <stdlib.h>
#include <cmath>
#include <iostream>

#include "rprop.h"
#include "gp_utils.h"

namespace libgp {

void RProp::init(double eps_stop, double Delta0, double Deltamin, double Deltamax, double etaminus, double etaplus) 
{
  this->Delta0   = Delta0;
  this->Deltamin = Deltamin;
  this->Deltamax = Deltamax;
  this->etaminus = etaminus;
  this->etaplus  = etaplus;
  this->eps_stop = eps_stop;

}

void RProp::maximize( GaussianProcess * gp, double RelTol, size_t n, bool verbose )
{
  int param_dim = gp->covf().get_param_dim();
  Eigen::VectorXd Delta = Eigen::VectorXd::Ones(param_dim) * Delta0;
  Eigen::VectorXd grad_old = Eigen::VectorXd::Zero(param_dim);
  Eigen::VectorXd params = gp->covf().get_loghyper();
  Eigen::VectorXd best_params = params;
  double best = log(0);


  double diff = 100;
  int cnt = 1;

  while( diff > RelTol && cnt < n ){

    Eigen::VectorXd grad = -gp->log_likelihood_gradient();
    grad_old = grad_old.cwiseProduct(grad);
    for (int j=0; j<grad_old.size(); ++j) {
      if (grad_old(j) > 0) {
        Delta(j) = std::min(Delta(j)*etaplus, Deltamax);        
      } else if (grad_old(j) < 0) {
        Delta(j) = std::max(Delta(j)*etaminus, Deltamin);
        grad(j) = 0;
      } 
      params(j) += -Utils::sign(grad(j)) * Delta(j);
    }
    grad_old = grad;
    if (grad_old.norm() < eps_stop) break;
    gp->covf().set_loghyper(params);
    double lik = gp->log_likelihood();
    
    
    if (lik > best) {
      best = lik;
      diff = (best_params-params).norm();
      best_params = params;
    }

    if (verbose) 
      printf("%d)  %.10le    (%.10le)\n ", cnt, -lik, diff);

    cnt ++;
  }
  
  gp->covf().set_loghyper(best_params);

}




}
