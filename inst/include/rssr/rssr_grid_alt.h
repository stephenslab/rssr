#ifndef RSSR_GRID_ALT_H
#define RSSR_GRID_ALT_H
#include <RcppEigen.h>
#include "rssr_types.h"

void rss_varbvsr_iter_alt(const c_Matrix_internal SiRiS,
                          const double sigma_beta,
                          const double logodds,
                          const c_arrayxd_internal betahat,
                          const c_arrayxd_internal se,
                          arrayxd_internal alpha,
                          arrayxd_internal mu,
                          arrayxd_internal SiRiSr,
                          bool reverse);


Rcpp::DataFrame grid_rss_varbvsr_alt(
    const c_Matrix_internal SiRiS,
    const c_arrayxd_internal sigma_beta,
    const c_arrayxd_internal logodds,
    const c_arrayxd_internal  betahat,
    const c_arrayxd_internal se,
    const c_arrayxd_internal talpha0,
    const c_arrayxd_internal tmu0,
    const c_arrayxd_internal tSiRiSr0,
    double tolerance,
    int itermax,
    bool verbose,
    bool lnz_tol);

Rcpp::DataFrame grid_search_rss_varbvsr_alt(
    const Matrix_external SiRis,
    const arrayxd_external sigma_beta,
    const arrayxd_external logodds,
    const arrayxd_external  betahat,
    const arrayxd_external  se,
    const arrayxd_external talpha0,
    const arrayxd_external tmu0,
    const arrayxd_external tSiRiSr0,
    double tolerance,
    int itermax,
    Rcpp::LogicalVector verbose,
    Rcpp::LogicalVector lnz_tol);




#endif
