#ifndef INCLUDE_BVSR
#define INCLUDE_BVSR
#include <RcppEigen.h>



void rss_varbvsr_update (const double betahat, const double se, const double sigma_beta, const Eigen::ArrayXd &SiRiS_snp, 
                         Eigen::ArrayXd &SiRiSr, const double SiRiSr_snp, const double logodds, double &alpha, double &mu);
void rss_varbvsr_iter(const Eigen::MappedSparseMatrix<double> SiRiS,
                      const double sigma_beta,
                      const double logodds,
                      const Eigen::Map<Eigen::ArrayXd> betahat,
                      const Eigen::Map<Eigen::ArrayXd> se,
                      Eigen::ArrayXd &alpha,
                      Eigen::ArrayXd &mu,
                      Eigen::ArrayXd &SiRiSr,
                      bool reverse);

#endif
