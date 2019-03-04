#include <fstream>
#include <stdio.h>

#include "gp_data.h"


gp_data::gp_data( string filename ){

  file.open( filename.c_str() );

  load_matrix( &file, data );

  Ns  = data.size();
  dim = data[0].size()-1;

}




gp_data::~gp_data( ){

}




void gp_data::split_train_test( int Ntr, int Ntst ){

  Ntrain = Ntr;
  Ntest  = Ntst;

  if( Ntrain+Ntest > Ns ){
  	cout << "Ntrain+Ntest must be less or equal than Ns (" << Ns << ")" << endl;
  	exit(1); // XXX better way?
  }

  // shuffle data and assign them to train and test vectors
  if(ind!=NULL) free(ind);
  ind = Utils::randperm(Ns);

}





void gp_data::set_gp( string covstr ){

  cov_str = covstr;

  delete gp;
  gp = new GaussianProcess( dim, covstr.c_str() );

  // initialize hyper parameter vector
  delete hprms;
  hprms = new Eigen::VectorXd( gp->covf().get_param_dim() );

  hprms->setZero();  
  *hprms = Eigen::VectorXd::Constant( gp->covf().get_param_dim(), 0.5 );

  // optimize hyperparameters
  gp->covf().set_loghyper( *hprms );


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

  double Ypred, tss=0., var;

  // total squared error
  for(int i = 0; i < Ntest; ++i) {

    double *Xtest = &data[ind[Ntrain+i]][0];

    Ypred = gp->f(Xtest);
    var   = gp->var(Xtest);

    for( int j=0; j<dim; j++) fp << Xtest[j] << "  ";

    fp << Ypred << "  " << var << "  " << data[ind[Ntrain+i]][dim] << endl;

    error = Ypred - data[ind[Ntrain+i]][dim] ;
    tss += error*error;
  }
  error = tss/Ntrain;


  // include also train data in the file
  for(int i = 0; i < Ntrain; ++i) {

    double *Xtrain = &data[ind[i]][0];

    Ypred = gp->f(Xtrain);
    var   = gp->var(Xtrain);

    for( int j=0; j<dim; j++) fp << Xtrain[j] << "  ";

    fp << Ypred << "  " << var << "  " << data[ind[i]][dim] << endl;

  }

  fp.close();

  return error;
}





double gp_data::get_error(){
  return error;
}




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
