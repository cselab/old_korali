// libgp - Gaussian process library for Machine Learning
// Copyright (c) 2013, Manuel Blum <mblum@informatik.uni-freiburg.de>
// All rights reserved.

#include <fstream>
#include <cmath>

#include "gp_wrap.h"

using namespace std;
using namespace Eigen;
using namespace libgp;
using namespace libgpwrap;

VectorXd test_gp( string datafile, string cov_str, string ident, grid gr, int Ntrain, int Ntest );

int main (int argc, char const *argv[]){

  static vector<string> cov_str = {
                                       "CovSEiso",
                                       "CovMatern3iso",
                                       "CovMatern5iso",
                                       "CovRQiso",
                                       "CovSEard",
                                       "CovSum ( CovMatern3iso, CovSEiso)",
                                       "CovSum ( CovSEiso, CovLinearard)",
                                       "CovSum ( CovSEiso, CovLinearone)",
                                       "CovProd( CovProd( CovLinearone, CovLinearone ), CovLinearone )"
                                  };

  // XXX InputDimFilter does not work!
  /*
  static vector<string> cov_str = {
                                  "CovSum( InputDimFilter(0/CovSEiso), InputDimFilter(1/CovSEiso) )"
                                  // " CovSum(CovSum(InputDimFilter(0/CovSEiso),InputDimFilter(1/CovSEiso)),InputDimFilter(2/CovSEiso))"
                                  };
  */

  srand48( 123 );
  
  Vector2d ab; ab << 0, 7;
  RowVectorXd xo;
  grid gr(ab,xo,1000,1);
  for(int i=0; i<cov_str.size(); i++){
    test_gp( "data/sincos.dat", cov_str[i], to_string(i+1), gr, 20, 100);
  }


  /*
  cout << "INFO: this example runs for long.. \n\n";
  Vector2d ab; ab << 0, 5;
  RowVectorXd xo(2); xo << 2.5, 2.5;
  grid gr(ab,xo,100,3);
  for(int i=0; i<cov_str.size(); i++){
      test_gp( "data/sincos3d.dat", cov_str[i], to_string(i+1), gr, 500, 100);
  }
  */

  return EXIT_SUCCESS;

}


VectorXd test_gp( string datafile, string cov_str, string ident, grid gr, int Ntrain, int Ntest ){

  static vector<string> filenames = { "data/gp", "data/out", "data/grad" };

  int verbose = 0;

  string str;

  gp_data gpd( datafile );

  gpd.split_train_test( Ntrain, Ntest, "random" );

  gpd.set_gp( cov_str, 0.5 );
  gpd.train_gp( 1e4, 1e-4, verbose );
  
  str = filenames[0] + ident + ".dat";
  gpd.get_gp().write( str.c_str() );

  str = filenames[1] + ident + ".dat";
  double gp_err = gpd.validate_gp( str );

  gpd.add_x( gr.X );
  gpd.eval_f();
  gpd.eval_df();

  str = filenames[2] + ident + ".dat";
  double gr_err = gpd.validate_df( gr.df_ind,  gr.h, str.c_str() );

  cout << "==========================================================\n";
  cout << cov_str << "\n\n";
  cout << "gp error   = " << gp_err << endl;
  cout << "grad error = " << gr_err << endl;
  cout << "gp parameters = " <<  gpd.get_gp().covf().get_loghyper().array().exp().transpose()  << endl;
  cout << "\n\n";

  VectorXd res(2);
  res << gp_err, gr_err;
  res *= sqrt(gr.h);
  return res;

}
