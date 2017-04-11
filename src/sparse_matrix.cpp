#include <RcppEigen.h>
#include "rssr.h"


// [[Rcpp::export]]
Eigen::SparseMatrix<double> SiRSi(const Eigen::MappedSparseMatrix<double> &R, const Eigen::VectorXd Si) {
  Eigen::SparseMatrix<double> m = Si.asDiagonal() *R*Si.asDiagonal();
  return(m);
}

 
// [[Rcpp::export]]
Eigen::MatrixXd SiRSi_d(const Matrix_external R, const Eigen::VectorXd Si) {
  return(Si.asDiagonal() *R*Si.asDiagonal());
}

// [[Rcpp::export]]
Eigen::SparseMatrix<double> genSymm(const Eigen::SparseMatrix<double> &R){
  Eigen::SparseMatrix<double> spId(R.rows(),R.cols());
  spId.reserve(R.rows());
  spId.setIdentity();
  Eigen::SparseMatrix<double> Rt = R.adjoint();
  Eigen::SparseMatrix<double> retmat=R+Rt-spId;
  return(retmat);
}
