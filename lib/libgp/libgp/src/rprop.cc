// libgp - Gaussian process library for Machine Learning
// Copyright (c) 2013, Manuel Blum <mblum@informatik.uni-freiburg.de>
// All rights reserved.

#include <stdlib.h>
#include <cmath>
#include <iostream>

#include "rprop.h"
#include "gp_utils.h"

namespace libgp {

void RProp::init(double dif_tol, double eps_stop, double Delta0, double Deltamin, double Deltamax, double etaminus, double etaplus) 
{
  this->Delta0   = Delta0;
  this->Deltamin = Deltamin;
  this->Deltamax = Deltamax;
  this->etaminus = etaminus;
  this->etaplus  = etaplus;
  this->eps_stop = eps_stop; //XXX sometimes the norm of the derivative is zero
  this->dif_tol  = dif_tol;

}

void RProp::maximize( GaussianProcess * gp, size_t n, bool verbose )
{
  int param_dim = gp->covf().get_param_dim();
  Eigen::VectorXd Delta = Eigen::VectorXd::Ones(param_dim) * Delta0;
  Eigen::VectorXd grad_old = Eigen::VectorXd::Zero(param_dim);
  Eigen::VectorXd params = gp->covf().get_loghyper();
  Eigen::VectorXd best_params = params;
  double best = log(0);


  int cnt = 1;
  double diff = 100;

  while( cnt <= n ){

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
    if (grad_old.norm() < eps_stop){
      printf("\nStop due to norm of gradient %.10le  <  %le\n\n",grad.norm(),eps_stop);
      break;
    }

    gp->covf().set_loghyper(params);
    double lik = gp->log_likelihood();
    
    
    if (lik > best) {
      best = lik;
      diff = (best_params-params).norm();
      best_params = params;
    }

   if( diff < dif_tol ){ 
      printf("\nStop due to norm of difference in parameters %.10le  <  %le\n\n", diff, dif_tol );
      break;
   }
    if (verbose) 
      printf("%d)  %.5le    (%.5le, %.5le)\n ", cnt, -lik, diff, grad.norm());

    cnt ++;
  }
  
  gp->covf().set_loghyper(best_params);

}




}
