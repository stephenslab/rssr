#ifndef RSSR_TBB_H
#define RSSR_TBB_H

#include <RcppEigen.h>
#include "rssr_types.h"

Rcpp::DataFrame grid_search_rss_varbvsr_tls(
    const  Matrix_external SiRiS,
    const arrayxd_external sigma_beta,
    const arrayxd_external logodds,
    const arrayxd_external  betahat,
    const arrayxd_external  se,
    const arrayxd_external talpha0,
    const arrayxd_external tmu0,
    const arrayxd_external tSiRiSr0,
    const double tolerance,
    const int itermax,
    Rcpp::LogicalVector lnz_tol,
    const int n,const int grainsize);

#endif
