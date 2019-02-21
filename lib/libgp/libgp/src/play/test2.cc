#include <Eigen/Dense>
int changeTwo(Eigen::Ref<Eigen::ArrayXd> f) {
for (int i = 0; i < 10; i++) {
f(i) = 2;
    }
return 0;
}
int main() {
    Eigen::MatrixXd randomMat = Eigen::MatrixXd::Random(1000,2);
changeTwo(randomMat.col(0));
}
