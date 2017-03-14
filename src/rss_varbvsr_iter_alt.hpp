#ifndef INCLUDE_BVSR_ALT
#define INCLUDE_BVSR_ALT
#include <RcppEigen.h>

void rss_varbvsr_iter_alt(const Eigen::SparseMatrix<double> &SiRiS,
                          const double sigma_beta,
                          const double logodds,
                          const Eigen::Map<Eigen::ArrayXd> betahat,
                          const Eigen::Map<Eigen::ArrayXd> se,
                          Eigen::ArrayXd &alpha,
                          Eigen::ArrayXd &mu,
                          Eigen::ArrayXd &SiRiSr,
                          bool reverse);

void compute_mu(const double betahat,
                const double se_square,
                const double sigma_square,
                const double alpha,
                double &mu,
                const double SiRiSr_snp);



void compute_alpha(const double sigma_square,
                   const double sigma_beta,
                   const double logodds,
                   const double mu,
                   double &alpha);

void compute_SiRiSr(const Eigen::SparseVector<double> &SiRiS_snp,
                    const  double r0,
                    const double r_new,
                    Eigen::VectorXd &SiRiSr);

Rcpp::List rss_varbvsr_naive_alt (const Eigen::MappedSparseMatrix<double> &SiRiS,
                                  const double sigma_beta,
                                  const double logodds,
                                  const Eigen::Map<Eigen::ArrayXd> betahat,
                                  const Eigen::Map<Eigen::ArrayXd> se,
                                  const Eigen::ArrayXd &talpha0,
                                  const Eigen::ArrayXd &tmu0,
                                  const Eigen::ArrayXd &tSiRiSr0,
                                  double tolerance,
                                  int itermax,
                                  Rcpp::LogicalVector verbose,
                                  Rcpp::LogicalVector lnz_tol);

#endif
