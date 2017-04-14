#ifndef KL_HPP
#define KL_HPP
#include <RcppEigen.h>
#include <limits>
#include "rssr_types.h"

typedef std::numeric_limits<double> double_lim;



inline Eigen::ArrayXd betavar(const c_arrayxd_internal p,const c_arrayxd_internal mu,const c_arrayxd_internal s);

double intklbeta_rssbvsr(const c_arrayxd_internal alpha,const c_arrayxd_internal mu,const c_arrayxd_internal sigma_square,double sigma_beta_square);

double intgamma(double logodds, const c_arrayxd_internal alpha);
double rel_err(double p0,double p1);
Eigen::ArrayXd rel_err(c_arrayxd_internal p0,c_arrayxd_internal p1);
double find_maxerr(const c_arrayxd_internal alpha,
                   const c_arrayxd_internal alpha0,
                   const c_arrayxd_internal r,
                   const c_arrayxd_internal r0);

double calculate_lnZ(const c_vectorxd_internal q,
                     const c_vectorxd_internal r,
                     const c_vectorxd_internal SiRiSr,
                     double logodds,
                     const c_vectorxd_internal sesquare,
                     const c_vectorxd_internal alpha,
                     const c_vectorxd_internal mu,
                     const c_vectorxd_internal s,
                     const double sigb);
double update_logodds(const c_arrayxd_internal alpha);

#endif
