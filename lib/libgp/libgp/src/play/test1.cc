#include <stdlib.h>
#include<iostream>
#include <Eigen/Dense>


using namespace std;
using namespace Eigen;



void sliceit( VectorXd     &x );
void sliceit( Ref<ArrayXd>  x );


int main(){


  cout << "Eigen: test slicing" << endl;

  MatrixXd A(5,5);
  
  A = MatrixXd::Identity(5,5);

  cout << A << "\n\n";

  cout << "--------------------\n";

  cout << A.col(0) << "\n\n";

  sliceit( A.col(0) );
  cout << endl;
  
  cout << A.col(0) << "\n\n";
  
  cout << "--------------------\n";

  cout << A << "\n\n";

  cout << "--------------------\n";


  VectorXd x = VectorXd::Zero(A.rows());
  sliceit(x);
  
  A.row(0) = x;

  cout << "--------------------\n";

  cout << A << "\n\n";

}






void sliceit( VectorXd &x ){

  x = VectorXd::Constant(x.size(),4);
  cout << "In sliceit:" << endl;
  cout << x << endl;

  return;
}




void sliceit(  Ref<ArrayXd> x ){

  x = ArrayXd::Constant(x.size(),2);

  cout << "In sliceit:" << endl;
  cout << x << endl;

  return;
}
