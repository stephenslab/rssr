#ifndef INCLUDE_BVSR_ALT
#define INCLUDE_BVSR_ALT
#include <RcppEigen.h>

void rss_varbvsr_iter_alt(const Eigen::SparseMatrix<double> SiRiS,
                          const double sigma_beta,
                          const double logodds,
                          const Eigen::ArrayXd betahat,
                          const Eigen::ArrayXd se,
                          Eigen::ArrayXd &alpha,
                          Eigen::ArrayXd &mu,
                          Eigen::ArrayXd &SiRiSr,
                          bool reverse);

#endif
