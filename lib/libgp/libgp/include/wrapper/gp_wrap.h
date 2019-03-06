#ifndef __GP_WRAP_H__
#define __GP_WRAP_H__

#include "gp_utils.h"
#include "rprop.h"
#include "cg.h"

#include <Eigen/Dense>
#include <string>
#include <vector>
#include <fstream>


namespace libgpwrap{

  using namespace libgp;
  using namespace std;
  using namespace Eigen;


  typedef std::vector< std::vector<double> > dmatrix;
  typedef std::vector<double> dvector;


  class grid{

    public:

      grid( Vector2d ab, RowVectorXd xo, int N, int ind );

      virtual ~grid(){};

      void display( );

      MatrixXd X;
      double h ;
      int df_ind ;
  };





  class gp_data {

    public:

      gp_data( string filename );

      virtual ~gp_data();

      void split_train_test( int Ntr, int Ntst, string tp="random");

      void set_gp( string covstr, double d=0 );

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

      double test_error = NAN;

      void load_matrix( istream* is, dmatrix & matrix );
  };


} // end of namespace libgpwrap



#endif // __GP_WRAP_H__
