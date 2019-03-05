#include <fstream>
#include <stdio.h>

#include "gp_wrap.h"


namespace libgpwrap{


  gp_data::gp_data( string filename ){

    file.open( filename.c_str() );

    load_matrix( &file, data );

    Ns  = data.size();
    dim = data[0].size()-1;

  }



  gp_data::~gp_data( ){
    delete gp;
    delete ind;
  }




  void gp_data::split_train_test( int Ntr, int Ntst, string tp ){

    Ntrain = Ntr;
    Ntest  = Ntst;

    if( Ntrain+Ntest > Ns ){
    	cout << "Ntrain+Ntest must be less or equal than Ns (" << Ns << ")" << endl;
    	exit(1); // XXX better way?
    }

    if(ind!=NULL) free(ind);

    if( !tp.compare("random") ){
      ind = Utils::randperm(Ns);
    }
    else{
      ind = new int[Ns];
      for(int i=0; i<Ns; i++) ind[i] = i;
    }

  }





  void gp_data::set_gp( string covstr, double d ){

    cov_str = covstr;

    delete gp;
    gp = new GaussianProcess( dim, covstr.c_str() );

    // initialize log-hyperparameter vector
    Eigen::VectorXd hprms = Eigen::VectorXd::Constant( gp->covf().get_param_dim(), d );
    gp->covf().set_loghyper( hprms );

    // add training patterns
    for(int i = 0; i < Ntrain; ++i){
      double *Xtrain = &data[ind[i]][0];
      double  Ytr = data[ind[i]][dim];
      gp->add_pattern( Xtrain, Ytr );
    }

  }




  void gp_data::train_gp( size_t nmax, double tol, int verbose){

    RProp rprop;
    rprop.init(tol);
    rprop.maximize( gp, nmax, verbose);

    //CG cg;
    //cg.maximize(&gp, 1000, 1);
  }




  GaussianProcess & gp_data::get_gp(){
    return *gp;
  }



  double gp_data::validate_gp( string filename ){

    ofstream fp;
    fp.open( filename.c_str() );

    fp << "%% X_1 | ... | X_N | prediction | variance | exact " << endl;

    double Ypred, var, tmp;

    // write test set and prediction in "filename"
    // and compute discrete normalized L2 norm on the test set
    test_error = 0;
    for( int i = 0; i < Ntest; ++i ){

      double *Xtest = &data[ind[Ntrain+i]][0];

      Ypred = gp->f(Xtest);
      var   = gp->var(Xtest);

      for( int j=0; j<dim; j++) fp << Xtest[j] << "  ";
      fp << Ypred << "  " << var << "  " << data[ind[Ntrain+i]][dim] << endl;

      tmp = Ypred - data[ind[Ntrain+i]][dim];
      test_error += tmp*tmp;
    }
    test_error = test_error/Ntrain;


    //  write also train set and prediction in "filename"
    for(int i = 0; i < Ntrain; ++i) {

      double *Xtrain = &data[ind[i]][0];

      Ypred = gp->f(Xtrain);
      var   = gp->var(Xtrain);

      for( int j=0; j<dim; j++) fp << Xtrain[j] << "  ";
      fp << Ypred << "  " << var << "  " << data[ind[i]][dim] << endl;
    }

    fp.close();

    return test_error;
  }




  double gp_data::get_error(){
    return test_error;
  }



  void gp_data::add_x( MatrixXd input ){

    int n = input.rows();
    int m = input.cols();

    assert(m==dim);

    X.resize(n,m);
    X << input;

  }


  void gp_data::eval_f( ){

    assert( X.rows()>0 && X.cols()==dim );

    F.resize( X.rows() );

    for( int i=0; i<X.rows(); i++){

      double x[dim];  // ugly; find better way. overload?
      for(int j=0 ; j<X.cols() ; j++) x[j] = X(i,j);

      F(i) = gp->f(x);
    }

  }



  void gp_data::eval_df( ){

    assert( X.rows()>0 && X.cols()==dim );

    DF.resize( X.rows(), X.cols() );

    for( int i=0; i<X.rows(); i++){

      double x[dim];  // ugly; find better way. overload?
      for(int j=0 ; j<X.cols() ; j++) x[j] = X(i,j);

      DF.row(i) = gp->dfdx( x );
    }

  }




  double gp_data::validate_df( int k, double h, string filename ){

    assert( X.rows()>0 && X.cols()==dim );
    assert( k < dim);

    int N = X.rows();

    DFFD.resize( X.rows()-2 );

    DFFD = ( F.tail(N-2) - F.head(N-2) )/(2*h);

    FILE *fp = fopen(filename.c_str(),"w");
    for( int i=1; i<N-1; i++ )
      fprintf(fp,"%.10le %.10le %.10le \n ",X(i,k),DF(i,k),DFFD(i-1));
    fclose(fp);

    return ( DFFD - DF.col(k).segment(1,N-2) ).norm();

  }







  //==============================================================================
  //==============================================================================

  // load matrix from an ascii text file.
  void gp_data::load_matrix( istream* is, dmatrix & matrix ){

      using namespace std;
  		const string& delim=" \t";

      string      line;
      string      strnum;

      // clear first
      matrix.clear();

      // parse line by line
      while (getline(*is, line))
      {
          matrix.push_back(vector<double>());

          for (string::const_iterator i = line.begin(); i != line.end(); ++ i)
          {
              // If i is not a delim, then append it to strnum
              if (delim.find(*i) == string::npos)
              {
                  strnum += *i;
                  if (i + 1 != line.end()) // If it's the last char, do not continue
                      continue;
              }

              // if strnum is still empty, it means the previous char is also a
              // delim (several delims appear together). Ignore this char.
              if (strnum.empty())
                  continue;

              // If we reach here, we got a number. Convert it to double.
              double number;

              istringstream(strnum) >> number;
              matrix.back().push_back(number);

              strnum.clear();
          }
      }
  }







  //==============================================================================
  //==============================================================================
  grid::grid( Vector2d ab, RowVectorXd xo, int N, int ind ){

    int dim = xo.size()+1;
    df_ind = ind-1;

    assert(df_ind<dim);

    h = ( ab(1)-ab(0) ) / (N-1);

    X.setZero(N,dim);
    X.leftCols(dim-1) = xo.replicate(N,1);
    X.col(dim-1) = VectorXd::LinSpaced( N, ab(0), ab(1));

    X.col(dim-1).swap(X.col(df_ind));
  }

  void grid::display(){
    cout << "\n\n" << X << "\n\n";
  }





}// end of namespace libgpwrap
