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


void rss_varbvsr_update(const double betahat,
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


void rss_varbvsr_iter(const c_Matrix_internal SiRiS,
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
  
  
  
  // Run coordinate ascent updates.
  // Repeat for each coordinate ascent update.
  size_t i=0;
  for (size_t j = 0; j < p; j++) {
    if(reverse){
      i=p-1-j;
    }else{
      i=j;
    }
    // Copy the kth column of matrix inv(S)*R*inv(S).
    // copySparseColumn(SiRiS_snp, k, SiRiS.elems, Ir, Jc, p);
    
    // Copy the kth element of vector inv(S)*R*inv(S)*r.
    double SiRiSr_snp = SiRiSr.coeff(i);
    
    // Perform the mean-field variational update.
    rss_varbvsr_update(betahat.coeff(i), se.coeff(i), sigma_beta, 
                       SiRiS.col(i), SiRiSr, SiRiSr_snp, 
                       logodds, alpha.coeffRef(i), mu.coeffRef(i));
  }
}



void rss_varbvsr_iter(const c_Matrix_internal SiRiS,
                      const c_arrayxd_internal sigma_beta,
                      const c_arrayxd_internal logodds,
                      const c_arrayxd_internal betahat,
                      const c_arrayxd_internal se,
                      arrayxxd_internal alpha,
                      arrayxxd_internal mu,
                      arrayxxd_internal SiRiSr,
                      bool reverse){
  
  
  // Get the number of SNPs (p) and coordinate ascent updates (m).
  const size_t p = betahat.size();
  const size_t gsize= sigma_beta.size();
  
  // Initialize outputs.
  
  // Store a single column of matrix inv(S)*R*inv(S).
  
  auto se_square = se*se;
  auto sigma_beta_square = sigma_beta*sigma_beta;
  Eigen::ArrayXXd sigma_square(p,gsize);
  for(size_t i=0; i<gsize; i++){
    sigma_square.col(i)= (se_square*sigma_beta_square.coeff(i))/(se_square+sigma_beta_square.coeff(i));
  }
  // Run coordinate ascent updates.
  // Repeat for each coordinate ascent update.
  size_t i=0;
  Eigen::ArrayXd r=alpha.row(0)*mu.row(0);
  Eigen::ArrayXd r_new=alpha.row(0)*mu.row(0);

  

  
  for (size_t j = 0; j < p; j++) {
    if(reverse){
      i=p-1-j;
    }else{
      i=j;
    }
    // Copy the kth column of matrix inv(S)*R*inv(S).
    // copySparseColumn(SiRiS_snp, k, SiRiS.elems, Ir, Jc, p);
    
    // Copy the kth element of vector inv(S)*R*inv(S)*r.
//    auto SiRiSr_snp = SiRiSr.row(i);
    
    // Compute the variational estimate of the posterior variance.
//    auto sigma_square = (se_square * sigma_beta_square) / (se_square + sigma_beta_square);
    
    // Update the variational estimate of the posterior mean.
    r = alpha.row(i) * mu.row(i);
    
    for(size_t k=0;k<gsize;k++){
       mu.coeffRef(i,k) = sigma_square.coeff(i,k) * (betahat.coeff(i) / se_square.coeff(i) + r.coeff(k) / se_square.coeff(i) - SiRiSr.coeff(i,k));
    }
    
    // Update the variational estimate of the posterior inclusion probability.
//    auto SSR = mu * mu / sigma_square;
    for(size_t k=0;k<gsize;k++){
      alpha.coeffRef(i,k) = sigmoid(logodds.coeff(k) + 0.5 * (log(sigma_square.coeff(i,k)/(sigma_beta_square.coeff(k))) + (mu.coeff(i,k)*mu.coeff(i,k))/sigma_square.coeff(i,k)));      
    }

    
    // Update SiRiSr = inv(S)*R*inv(S)*r
    r_new = alpha.row(i) * mu.row(i);
    //  cblas_daxpy(p,(r_new-r),SiRiS_snp.data(),1,SiRiSr.data(),1);
    for(int k=0; k<gsize;k++){
      SiRiSr.col(k)=(SiRiS.col(k)).array()*(r_new.coeff(k)-r.coeff(k));    
    }
  }
}






