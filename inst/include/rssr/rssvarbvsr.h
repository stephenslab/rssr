#ifndef INCLUDE_BVSR
#define INCLUDE_BVSR
#include <RcppEigen.h>
#include "rssr_types.h"
using namespace Eigen;


void rss_varbvsr_update(const double betahat,
                        const double se_square,
                        const double sigma_beta_square,
                        const c_arrayxd_internal SiRiS_snp,
                        const double sigma_square,
                        arrayxd_internal SiRiSr,
                        const double SiRiSr_snp,
                        const double ssrat,
                        const double logodds,
                        double &alpha,
                        double &mu);

void rss_varbvsr_iter(const c_sparseMatrix_internal SiRiS,
                         const double sigma_beta_square,
                         const c_arrayxd_internal sigma_square,
                         const double logodds,
                         const c_arrayxd_internal betahat,
                         const c_arrayxd_internal se_square,
                         const c_arrayxd_internal ssrat,
                         arrayxd_internal alpha,
                         arrayxd_internal mu,
                         arrayxd_internal SiRiSr,
                         bool reverse);

void rss_varbvsr_iter(const c_Matrix_internal SiRiS,
                      const double sigma_beta_square,
                      const c_arrayxd_internal sigma_square,
                      const double logodds,
                      const c_arrayxd_internal betahat,
                      const c_arrayxd_internal se_square,
                      const c_arrayxd_internal ssrat,
                      arrayxd_internal alpha,
                      arrayxd_internal mu,
                      arrayxd_internal SiRiSr,
                      bool reverse);

// void rss_varbvsr_iter(const c_Matrix_internal SiRiS,
//                       const c_arrayxd_internal sigma_beta_square,
//                       const c_arrayxxd_internal sigma_square,
//                       const c_arrayxd_internal logodds,
//                       const c_arrayxd_internal betahat,
//                       const c_arrayxd_internal se_square,
//                       const c_arrayxxd_internal ssrat,
//                       arrayxxd_internal alpha,
//                       arrayxxd_internal mu,
//                       arrayxxd_internal SiRiSr,
//                       bool reverse);



#endif
