#include <RcppEigen.h>
#include "rssr.h"
//#include "mkl.h"
#include <math.h>

using namespace Rcpp;



// RcppExport SEXP start_profiler(SEXP str) {
//   ProfilerStart(as<const char*>(str));
//   return R_NilValue;
// }
// 
// RcppExport SEXP stop_profiler() {
//   ProfilerStop();
//   return R_NilValue;
// }

// USAGE: run a single coordinate ascent of variational update to fit RSS-BVSR


void rss_varbvsr_update (const double betahat,
                         const double se,
                         const double sigma_beta,
                         const c_arrayxd_internal SiRiS_snp,
                         arrayxd_internal  SiRiSr,
                         const double SiRiSr_snp,
                         const double logodds,
                         double &alpha,
                         double &mu) {
  
  double se_square = se * se;
  double sigma_beta_square = sigma_beta * sigma_beta;
  size_t p=SiRiS_snp.size();
  // Compute the variational estimate of the posterior variance.
  double sigma_square = (se_square * sigma_beta_square) / (se_square + sigma_beta_square);
  
  // Update the variational estimate of the posterior mean.
  double r = alpha * mu;
  mu = sigma_square * (betahat / se_square + r / se_square - SiRiSr_snp);
  
  // Update the variational estimate of the posterior inclusion probability.
  double SSR = mu * mu / sigma_square;
  alpha = sigmoid(logodds + 0.5 * (log(sigma_square/(sigma_beta_square)) + SSR));

  // Update SiRiSr = inv(S)*R*inv(S)*r
  double r_new = alpha * mu;
//  cblas_daxpy(p,(r_new-r),SiRiS_snp.data(),1,SiRiSr.data(),1);
  SiRiSr+=(SiRiS_snp*(r_new-r));
}













void rss_varbvsr_iter(const c_sparseMatrix_internal SiRiS,
                      const double sigma_beta,
                      const double logodds,
                      const c_arrayxd_internal betahat,
                      const c_arrayxd_internal se,
                      arrayxd_internal alpha,
                      arrayxd_internal mu,
                      arrayxd_internal SiRiSr,
                      bool reverse){
  
  
  // Get the number of SNPs (p) and coordinate ascent updates (m).
  const size_t p = betahat.size();
  
  // Initialize outputs.
  
  // Store a single column of matrix inv(S)*R*inv(S).
  
  Eigen::ArrayXd  SiRiS_snp(p);
  Eigen::VectorXd  SiRiS_snp_v(p);
  
  // Run coordinate ascent updates.
  // Repeat for each coordinate ascent update.
  size_t i=0;
  for (size_t j = 0; j < p; j++) {
    if(reverse){
      i=p-1-j;
    }else{
      i=j;
    }
    SiRiS_snp_v=(SiRiS.col(i));
    SiRiS_snp=SiRiS_snp_v.array();
    // Copy the kth column of matrix inv(S)*R*inv(S).
    // copySparseColumn(SiRiS_snp, k, SiRiS.elems, Ir, Jc, p);
    
    // Copy the kth element of vector inv(S)*R*inv(S)*r.
    double SiRiSr_snp = SiRiSr.coeff(i);
    
    // Perform the mean-field variational update.
    rss_varbvsr_update(betahat.coeff(i), se.coeff(i), sigma_beta, 
                       SiRiS_snp, SiRiSr, SiRiSr_snp, 
                       logodds, alpha.coeffRef(i), mu.coeffRef(i));
  }
}




