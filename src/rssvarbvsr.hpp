#ifndef INCLUDE_BVSR
#define INCLUDE_BVSR
#include <RcppEigen.h>





//typedef MSpMat::InnerIterator InIterMat;
//typedef Eigen::SparseVector<double> SpVec;
//typedef SpVec::InnerIterator InIterVec;


void rss_varbvsr_update (const double betahat, const double se, const double sigma_beta, const Eigen::ArrayXd &SiRiS_snp, 
                         Eigen::ArrayXd &SiRiSr, const double SiRiSr_snp, const double logodds, double &alpha, double &mu);
void rss_varbvsr_iter(const Eigen::SparseMatrix<double> SiRiS,
                      const double sigma_beta,
                      const double logodds,
                      const Eigen::ArrayXd betahat,
                      const Eigen::ArrayXd se,
                      Eigen::ArrayXd &alpha,
                      Eigen::ArrayXd &mu,
                      Eigen::ArrayXd &SiRiSr,
                      bool reverse);
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
