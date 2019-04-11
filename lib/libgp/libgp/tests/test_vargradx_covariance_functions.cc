// libgp - Gaussian process library for Machine Learning
// Copyright (c) 2011, Manuel Blum <mblum@informatik.uni-freiburg.de>
// All rights reserved.

#include "cov_factory.h"
#include "gp_wrap.h"

#include <Eigen/Dense>
#include <gtest/gtest.h>


using ::testing::TestWithParam;
using ::testing::Values;

using namespace libgpwrap;
using namespace Eigen;

class Vargrad_x_Test : public TestWithParam<std::string> {

  protected:

    gp_data *gpd;
    grid *gr;

    virtual void SetUp(){
      srand48( 123 );
      int Ntrain = 10, Ntest = 10;

      gpd = new gp_data( "data/sincos.dat" );
      gpd->split_train_test( Ntrain, Ntest, "random" );
      gpd->set_gp( GetParam() );
      gpd->train_gp( 1e4, 1e-4, 0 );

      Vector2d ab; ab << 0, 7;
      RowVectorXd xo;
      gr = new grid(ab,xo,1000,1);

      gpd->add_x( gr->X );
      gpd->eval_var();
      gpd->eval_dvar();

    }

    virtual void TearDown(){
      delete gpd;
      delete gr;
    }

    double vargradx_error( ){
      return gpd->validate_dvar( gr->df_ind,  gr->h, "" );
    }

};



TEST_P( Vargrad_x_Test, EqualToNumerical ){
  
    EXPECT_TRUE( fabs(vargradx_error()) < 5e-2 );

}



INSTANTIATE_TEST_SUITE_P( CovarianceFunction, Vargrad_x_Test, Values(
/*          "CovLinearard",
          "CovLinearone",
          "CovMatern3iso",
          "CovMatern5iso",
          "CovProd(CovSEiso, CovMatern3iso)",
          "CovRQiso",
          "CovSEard", */
//          "CovSEiso",
          "CovSum(CovSEiso, CovNoise)"
/*          "CovSum(CovLinearard, CovNoise)",
          "InputDimFilter(0/CovSEiso)",
          "InputDimFilter(0/CovSum(CovSEiso, CovNoise))" */
          ));



//==========================================================================================================================

class Vargrad_x_nd_Test : public TestWithParam<std::string> {

  protected:

    gp_data *gpd;
    grid *gr;
    RowVectorXd xo;

    virtual void SetUp(){
      srand48( 123 );
      int Ntrain = 10, Ntest = 10;

      gpd = new gp_data( "data/sincos3d.dat" );
      gpd->split_train_test( Ntrain, Ntest, "random" );
      gpd->set_gp( GetParam() );
      gpd->train_gp( 1e4, 1e-4, 0 );
    }

    virtual void TearDown(){
      delete gpd;
      delete gr;
    }

    double vargradx_error( int k ){
      Vector2d ab; ab << 0, 5;
      RowVectorXd xo(2); xo << 1.5, 1.5;
      gr = new grid(ab,xo,1000,k);
      gpd->add_x( gr->X );
      gpd->eval_var();
      gpd->eval_dvar();
      return gpd->validate_dvar( gr->df_ind,  gr->h, "dvar_test.out" );
    }

};



TEST_P( Vargrad_x_nd_Test, EqualToNumerical ){
  for(int k=1; k<=3; k++){
    //EXPECT_TRUE( fabs(vargradx_error(k)) < 0.005 );
  }
}



INSTANTIATE_TEST_SUITE_P( CovarianceFunction, Vargrad_x_nd_Test, Values(
          "CovLinearard",
          "CovLinearone",
          "CovMatern3iso",
          "CovMatern5iso",
          "CovProd(CovSEiso, CovMatern3iso)",
          "CovRQiso",
          "CovSEard",
          "CovSEiso",
          "CovSum(CovSEiso, CovNoise)",
          "CovSum(CovLinearard, CovNoise)",
          "InputDimFilter(0/CovSEiso)",
          "InputDimFilter(0/CovSum(CovSEiso, CovNoise))"
));
