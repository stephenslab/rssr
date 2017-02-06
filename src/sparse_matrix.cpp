#include <RcppEigen.h>



// [[Rcpp::export]]
Eigen::SparseMatrix<double> SiRSi(const Eigen::SparseMatrix<double> &R, const Eigen::VectorXd Si) {
  Eigen::SparseMatrix<double> m = Si.asDiagonal() *R*Si.asDiagonal();
  return(m);
}
