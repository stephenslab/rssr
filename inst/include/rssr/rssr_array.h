#ifndef RSSR_ARRAY_H
#define RSSR_ARRAY_H

#include <RcppEigen.h>
#include "rssr_types.h"


Eigen::ArrayXd rss_varbvsr_squarem_iter_array(const c_Matrix_internal SiRiS,
                                              const c_arrayxd_internal sigma_beta,
                                              const c_arrayxd_internal logodds,
                                              const c_arrayxd_internal  betahat,
                                              const c_arrayxd_internal se,
                                              const c_arrayxd_internal talpha0,
                                              const c_arrayxd_internal tmu0,
                                              const c_arrayxd_internal tSiRiSr0,
                                              double tolerance,
                                              int itermax,
                                              Rcpp::LogicalVector lnz_tol);






#endif
