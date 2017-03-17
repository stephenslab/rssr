#ifndef RSSR_GRID_ALT_GRID_H
#define RSSR_GRID_ALT_GRID_H
#include <RcppEigen.h>
#include "rssr_types.h"

Rcpp::DataFrame rss_varbvsr_alt_grid(const c_sparseMatrix_internal SiRiS,
                                     const c_arrayxd_internal sigma_beta,
                                     const c_arrayxd_internal logodds,
                                     const c_arrayxd_internal betahat,
                                     const c_arrayxd_internal se,
                                     const c_arrayxd_internal talpha0,
                                     const c_arrayxd_internal tmu0,
                                     const c_arrayxd_internal tSiRiSr0,
                                     double tolerance,
                                     int itermax,
                                     bool verbose,
                                     bool lnztol);
#endif
