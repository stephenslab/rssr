#ifndef INCLUDE_BVSR_ALT
#define INCLUDE_BVSR_ALT
#include <RcppEigen.h>
#include "rssr_types.h"



Eigen::SparseMatrix<double> SiRSi(const Eigen::MappedSparseMatrix<double> &R, const Eigen::VectorXd Si);

Eigen::MatrixXd SiRSi_d(const Matrix_external R, const Eigen::VectorXd Si);

#endif
