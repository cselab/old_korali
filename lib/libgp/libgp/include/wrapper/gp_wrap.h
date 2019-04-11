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

  typedef std::vector< std::vector<double> > dmatrix;
  typedef std::vector<double> dvector;

  class grid{

    public:

      grid( Eigen::Vector2d ab, Eigen::RowVectorXd xo, int N, int ind );

      virtual ~grid(){};

      void display( );

      Eigen::MatrixXd X;
      double h ;
      int df_ind ;
  };


  class gp_data {

    public:

      gp_data( std::string filename );

      virtual ~gp_data();

      void split_train_test( int Ntr, int Ntst, std::string tp="random");

      void set_gp( std::string covstr, double d=0 );

      void train_gp( size_t Nmax, double tol, int verbose=1);

      libgp::GaussianProcess & get_gp();

      double validate_gp( std::string filename );

      double get_error();

      void add_x( Eigen::MatrixXd input );

      void eval_f( );

      void eval_df( );
      
      void eval_var( );
      
      void eval_dvar( );

      double validate_df( int k, double h, std::string filename );
      
      double validate_dvar( int k, double h, std::string filename );

    protected:

      std::ifstream file;

      dmatrix data;

      Eigen::MatrixXd X, DF, DV;
      
      Eigen::VectorXd F, DFFD, V, DVVD;

      int Ns, dim;

      int Ntrain, Ntest;

      int *ind = NULL;

      std::string cov_str;

      libgp::GaussianProcess *gp = nullptr;

      double test_error = NAN;

      void load_matrix( std::istream* is, dmatrix & matrix );
  };


} // end of namespace libgpwrap



#endif // __GP_WRAP_H__
