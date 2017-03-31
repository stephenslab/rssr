#include <RcppEigen.h>
#include "rssr.h"


//' Single update of RSS with variational method
//' This function is a very close translation of the original implementation of RSS. It is kept here for testing purposes
//' It performs a single update
// [[Rcpp::export]]
Eigen::MatrixXd rss_varbvsr_iter_naive_reference (Eigen::SparseMatrix<double> SiRiS,
                                                  Rcpp::NumericVector sigma_beta,
                                                  const Eigen::VectorXd logodds,
                                                  const Eigen::VectorXd betahat,
                                                  const Eigen::VectorXd se,
                                                  const Eigen::ArrayXd alpha0,
                                                  const Eigen::ArrayXd mu0,
                                                  const Eigen::ArrayXd SiRiSr0,bool reverse){
 
  
  // Get the number of SNPs (p) and coordinate ascent updates (m).
  const size_t p = betahat.size();
  
  // Initialize outputs.
  
  // Store a single column of matrix inv(S)*R*inv(S).
  
  Eigen::ArrayXd alpha=alpha0;
  Eigen::ArrayXd mu=mu0;
  Eigen::ArrayXd SiRiSr=SiRiSr0;
  
  Eigen::ArrayXd  SiRiS_snp(p);
  Eigen::VectorXd  SiRiS_snp_v(p);
  Eigen::MatrixXd retmat(p,3);
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
    double SiRiSr_snp = SiRiSr(i);
    
    // Perform the mean-field variational update.
    rss_varbvsr_update(betahat(i), se(i), sigma_beta[i], 
                       SiRiS_snp, SiRiSr, SiRiSr_snp, 
                       logodds(i), alpha(i), mu(i));
  }
  retmat << alpha,mu,SiRiSr;
  return(retmat);
  // Free the dynamically allocated memory.
}
