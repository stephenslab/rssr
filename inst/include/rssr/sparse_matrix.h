#ifndef INCLUDE_BVSR_ALT
#define INCLUDE_BVSR_ALT
#include <RcppEigen.h>



Eigen::SparseMatrix<double> SiRSi(const Eigen::MappedSparseMatrix<double> &R, const Eigen::VectorXd Si);

#endif
