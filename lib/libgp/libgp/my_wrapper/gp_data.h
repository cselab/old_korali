#include "gp_utils.h"
#include "rprop.h"
#include "cg.h"

#include <Eigen/Dense>
#include <string>
#include <vector>

#ifndef __GP_DATA_H__
#define __GP_DATA_H__

  using namespace libgp;
  using namespace std;
  using namespace Eigen;

  typedef std::vector< std::vector<double> > dmatrix;
  typedef std::vector<double> dvector;


  class gp_data {

    public:

      gp_data( string filename );

      virtual ~gp_data ();

      void split_train_test( int Ntr, int Ntst );

      void set_gp( string covstr );

      void train_gp( size_t Nmax, double tol, int verbose=1);

      GaussianProcess & get_gp();

      double validate_gp( string filename );

      double get_error();

      void add_x( MatrixXd input );

      void eval_f( );

      void eval_df( );

      double validate_df( int k, double h, string filename );



    protected:

      ifstream file;

      dmatrix data;

      MatrixXd X, DF;
      VectorXd F, DFFD;

      int Ns, dim;

      int Ntrain, Ntest;

      int *ind = NULL;

      string cov_str;

      GaussianProcess *gp = nullptr;

      Eigen::VectorXd *hprms = nullptr;

      double error = NAN;

      void load_matrix( istream* is, dmatrix & matrix );


  };


#endif // __GP_DATA_H__
