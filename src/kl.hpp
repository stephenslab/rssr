#ifndef KL_HPP
#define KL_HPP
#include <RcppEigen.h>
#include <limits>

typedef std::numeric_limits<double> double_lim;

Eigen::ArrayXd betavar(const Eigen::ArrayXd &p,const Eigen::ArrayXd &mu,Eigen::ArrayXd &s);

double intklbeta_rssbvsr(const Eigen::ArrayXd &alpha,const Eigen::ArrayXd &mu,const Eigen::ArrayXd &sigma_square,double sigma_beta_square);

double intgamma(double logodds, const Eigen::ArrayXd &alpha);
double rel_err(double p0,double p1);
double find_maxerr(const Eigen::ArrayXd &alpha,
                   const Eigen::ArrayXd &alpha0,
                   const Eigen::ArrayXd &r,
                   const Eigen::ArrayXd &r0);

double calculate_lnZ(const Eigen::VectorXd &q,
                     const Eigen::VectorXd &r,
                     const Eigen::VectorXd &SiRiSr,
                     double logodds,
                     const Eigen::VectorXd &sesquare,
                     const Eigen::VectorXd &alpha,
                     const Eigen::VectorXd &mu,
                     const Eigen::VectorXd &s,
                     const double sigb);

#endif