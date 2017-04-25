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
  double sigma_beta_square=sigma_beta*sigma_beta;
  size_t p=SiRiS_snp.size();
  // Compute the variational estimate of the posterior variance.
  double sigma_square = (se_square * sigma_beta_square) / (se_square + sigma_beta_square);
  double *SiRiS_p;
  // *SiRiS_p=SiRiS_snp.data();
  
  // Update the variational estimate of the posterior mean.
  double r = alpha * mu;
  // double tmu=mu;
  mu = sigma_square * (betahat / se_square + r / se_square - SiRiSr_snp);
  // if(!std::isfinite(mu)){
  //   Rcpp::Rcerr<<"old mu was "<<tmu<<std::endl;
  //   Rcpp::Rcerr<<"new mu is "<<mu<<std::endl;
  //   Rcpp::stop("mu is not finite");
  // }
  
  // Update the variational estimate of the posterior inclusion probability.
  double SSR = mu * mu / sigma_square;
  // double talpha=alpha;
  
  alpha = sigmoid(logodds + 0.5 * (log(sigma_square)-log(sigma_beta_square) + SSR));
  // if(!std::isfinite(alpha)){
  //   double psigmoid =logodds+0.5*(log(sigma_square)-log(sigma_beta_square)+SSR);
  //   Rcpp::Rcerr<<"old alpha was "<<talpha<<std::endl;
  //   Rcpp::Rcerr<<"new alpha is "<<alpha<<std::endl;
  //   Rcpp::Rcerr<<"psigmoid is "<<psigmoid<<std::endl;
  //   Rcpp::Rcerr<<"logodds is "<<logodds<<std::endl;
  //   Rcpp::Rcerr<<"sigma_square is "<<logodds<<std::endl;
  //   Rcpp::Rcerr<<"sigma_beta_square is "<<sigma_beta_square<<std::endl;
  //   Rcpp::Rcerr<<"sigma_beta is "<<sigma_beta<<std::endl;
  //   Rcpp::Rcerr<<"SSR is "<<SSR<<std::endl;
  //   
  //   
  //   Rcpp::stop("alpha is not finite");
  // }
  // Update SiRiSr = inv(S)*R*inv(S)*r
  double r_new = alpha * mu-r;
  cblas_daxpy(p,r_new,SiRiS_snp.data(),1,SiRiSr.data(),1);
  //  cblas_daxpy(p,(r_new-r),SiRiS_snp.data(),1,SiRiSr.data(),1);
  //  SiRiSr+=(SiRiS_snp*(r_new-r));
}



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
                        double &mu) {

  
  const size_t p=SiRiS_snp.size();

  
  double r = alpha * mu;
  
  mu = sigma_square * (betahat / se_square + r/se_square - SiRiSr_snp);
  alpha = sigmoid(logodds + 0.5 * (ssrat + mu*mu/sigma_square));
  
  double r_new = alpha * mu-r;
  cblas_daxpy(p,r_new,SiRiS_snp.data(),1,SiRiSr.data(),1);
}






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
                      bool reverse){
  
  mkl_set_num_threads_local(1);
  size_t p=betahat.size();
  // Get the number of SNPs (p) and coordinate ascent updates (m).
  
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

    
    // Perform the mean-field variational update.
    rss_varbvsr_update(betahat.coeff(i),
                         se_square.coeff(i),
                         sigma_beta_square, 
                         SiRiS.col(i),
                         sigma_square.coeff(i),
                         SiRiSr,
                         SiRiSr.coeff(i),
                         ssrat.coeff(i),
                         logodds,
                         alpha.coeffRef(i),
                         mu.coeffRef(i));
  }
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







