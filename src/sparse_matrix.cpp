#include <RcppEigen.h>
#include "rssr.h"


// [[Rcpp::export]]
Eigen::SparseMatrix<double> SiRSi(const Eigen::MappedSparseMatrix<double> &R, const Eigen::VectorXd Si) {
  if(R.rows()!=Si.size()){
    Rcpp::stop("Si must have as many elements as R has rows.");
  }
  Eigen::SparseMatrix<double> m = Si.asDiagonal() *R*Si.asDiagonal();
  return(m);
}

//[[Rcpp::export]]
Eigen::MatrixXd SiRSi_c(const Matrix_external R, const Eigen::ArrayXd se, const Eigen::ArrayXd betahat,const double n){
  if(R.rows()!=se.size()){
    Rcpp::stop("se must have as many elements as R has rows.");
  }
  Eigen::ArrayXd svec=(se.square()+betahat.square()/n).sqrt().inverse();
  return(svec.matrix().asDiagonal()*R*svec.matrix().asDiagonal());
}

// [[Rcpp::export]]
Eigen::MatrixXd SiRSi_a(const Matrix_external R, const Eigen::VectorXd Si) {
  if(R.rows()!=Si.size()){
    Rcpp::stop("Si must have as many elements as R has rows.");
  }
  return((R.array().colwise()*Si.array()).rowwise()*Si.transpose().array());
}

 
 // [[Rcpp::export]]
 Eigen::MatrixXd SiRSi_d(const Matrix_external R, const Eigen::VectorXd Si) {
   if(R.rows()!=Si.size()){
     Rcpp::stop("Si must have as many elements as R has rows.");
   }
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
