// libgp - Gaussian process library for Machine Learning
// Copyright (c) 2013, Manuel Blum <mblum@informatik.uni-freiburg.de>
// All rights reserved.

#include "gp_utils.h"
#include "rprop.h"
#include "cg.h"

#include <Eigen/Dense>
#include <fstream>
#include <string>
#include <vector>


using namespace libgp;
using namespace std;



void load_matrix( istream* is, vector< vector<double> >* matrix );




int main (int argc, char const *argv[])
{

  if (argc < 4){
  	cerr << "Usage: " << argv[0] << " Ntrain  Ntest  DataFile" << endl;
  	return 1;
  }


  int Ntrain = atoi(argv[1]);
  int Ntest  = atoi(argv[2]);


  // open file and load the data
  ifstream file;
  file.open( argv[3]);

  vector< vector<double> > data;
  load_matrix( &file, &data );

  int Ns  = data.size();
  int dim = data[0].size()-1;

  if( Ntrain+Ntest > Ns ){
  	cout << "Ntrain+Ntest must be less or equal than Ns (" << Ns << ")" << endl;
  	return 1;
  }


  // shuffle data and assign them to train and test vectors
  int *ind = Utils::randperm(Ns);

  vector< vector<double> >    xtrain( Ntrain, vector<double>(dim));
  vector< vector<double> >    xtest(  Ntest,  vector<double>(dim));

  vector<double>  ytrain(Ntrain);
  vector<double>  ytest( Ntest);

  for( int i=0; i<Ntrain; i++ ){
    for( int j=0; j<dim; j++ )
  		xtrain[i][j] = data[ind[i]][j];
  	ytrain[i]=data[ind[i]][dim];
  }

  for( int i=0; i<Ntest; i++ ){
  	for( int j=0; j<dim; j++ )
  		xtest[i][j] = data[ind[Ntrain+i]][j];
  	ytest[i]=data[ind[Ntrain+i]][dim];
  }





  // initialize Gaussian process for dim-D input using the squared exponential
  // covariance function with additive white noise.
  GaussianProcess gp( dim, "CovSum ( CovSEiso, CovNoise)" );

  // initialize hyper parameter vector
  Eigen::VectorXd   params( gp.covf().get_param_dim() );
  cout << "Rprop optimizer" << endl;
  params.setZero();

  // add training patterns
  double *Ytrain = &ytrain[0];
  for(int i = 0; i < Ntrain; ++i){
    double *Xtrain = &xtrain[i][0];
    gp.add_pattern( Xtrain, Ytrain[i] );
  }


  // optimize hyperparameters
  gp.covf().set_loghyper( params );

  RProp rprop;
  rprop.init(1e-8);
  rprop.maximize(&gp, 100, 1);


  // write data to file
  ofstream fp;
  fp.open ("out.txt");
  fp << "%% X_1 | ... | X_N | prediction | variance | exact " << endl;

  double Ypred, tss=0., error, var;

  // total squared error
  for(int i = 0; i < Ntest; ++i) {

    double *Xtest = &xtest[i][0];

    Ypred = gp.f(Xtest);
    var   = gp.var(Xtest);

    for( int j=0; j<dim; j++) fp << Xtest[j] << "  ";

    fp << Ypred << "  " << var << "  " << ytest[i] << endl;

    error = Ypred - ytest[i] ;
    tss += error*error;
  }

  // include also train data in the file
  for(int i = 0; i < Ntrain; ++i) {

    double *Xtrain = &xtrain[i][0];

    Ypred = gp.f(Xtrain);
    var   = gp.var(Xtrain);

    for( int j=0; j<dim; j++) fp << Xtrain[j] << "  ";

    fp << Ypred << "  " << var << "  " << ytrain[i] << endl;

    error = Ypred - ytrain[i] ;
  }







  fp.close();

  gp.write("gp.txt");

  cout << endl <<  "test set mse = " << tss/Ntest << endl;

  cout << endl << "kernel hyper-parameters:" << endl;
  cout << gp.covf().get_loghyper();

  cout << endl << endl;


  return EXIT_SUCCESS;

}















// load matrix from an ascii text file.
void load_matrix(istream* is,
				vector< vector<double> >* matrix)
{
    using namespace std;
		const string& delim=" \t";

    string      line;
    string      strnum;

    // clear first
    matrix->clear();

    // parse line by line
    while (getline(*is, line))
    {
        matrix->push_back(vector<double>());

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
            double       number;

            istringstream(strnum) >> number;
            matrix->back().push_back(number);

            strnum.clear();
        }
    }
}
