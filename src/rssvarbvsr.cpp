#include <RcppEigen.h>
#include "rssr.h"

#include <math.h>

using namespace Rcpp;


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
  SiRiSr+=SiRiS_snp*r_new;
  
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
  
  // mkl_set_num_threads_local(1);
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
    
    
    
    rss_varbvsr_update(betahat.coeff(i),
                       se_square.coeff(i),
                       sigma_beta_square, 
                       SiRiS_snp,
                       sigma_square.coeff(i),
                       SiRiSr,
                       SiRiSr.coeff(i),
                       ssrat.coeff(i),
                       logodds,
                       alpha.coeffRef(i),
                       mu.coeffRef(i));
  }
}




//[[Rcpp::export]]
Rcpp::List wrap_rss_varbvsr_iter(const Matrix_external SiRiS,
                                 const double sigma_beta,
                                 const double logodds,
                                 const arrayxd_external betahat,
                                 const arrayxd_external se,
                                 const arrayxd_external alpha,
                                 const arrayxd_external mu,
                                 const arrayxd_external SiRiSr,
                                 bool reverse){
  
  Eigen::ArrayXd talpha=alpha;
  Eigen::ArrayXd tmu=mu;
  Eigen::ArrayXd tSiRiSr=SiRiSr;
  double sigma_beta_square=sigma_beta*sigma_beta;
  Eigen::ArrayXd sesquare=se.square();
  Eigen::ArrayXd  s= (sesquare*(sigma_beta*sigma_beta))/(sesquare+(sigma_beta*sigma_beta));
  Eigen::ArrayXd ssrat((s/sigma_beta_square).log());
  rss_varbvsr_iter(SiRiS,sigma_beta_square,s,logodds,betahat,sesquare,ssrat,talpha,tmu,tSiRiSr,reverse);
  return Rcpp::List::create(Rcpp::Named("alpha1")=talpha,
                            Rcpp::Named("mu1")=tmu,
                            Rcpp::Named("SiRiSr")=tSiRiSr);
}



//[[Rcpp::export]]
Rcpp::List wrap_rss_varbvsr_iter_sp(const sparseMatrix_external SiRiS,
                                    const double sigma_beta,
                                    const double logodds,
                                    const arrayxd_external betahat,
                                    const arrayxd_external se,
                                    const arrayxd_external alpha,
                                    const arrayxd_external mu,
                                    const arrayxd_external SiRiSr,
                                    bool reverse){
  Eigen::ArrayXd talpha=alpha;
  Eigen::ArrayXd tmu=mu;
  Eigen::ArrayXd tSiRiSr=SiRiSr;
  double sigma_beta_square=sigma_beta*sigma_beta;
  Eigen::ArrayXd sesquare=se.square();
  Eigen::ArrayXd  s= (sesquare*(sigma_beta*sigma_beta))/(sesquare+(sigma_beta*sigma_beta));
  Eigen::ArrayXd ssrat((s/sigma_beta_square).log());
  rss_varbvsr_iter(SiRiS,sigma_beta_square,s,logodds,betahat,sesquare,ssrat,talpha,tmu,tSiRiSr,reverse);
  return Rcpp::List::create(Rcpp::Named("alpha1")=talpha,
                            Rcpp::Named("mu1")=tmu,
                            Rcpp::Named("SiRiSr")=tSiRiSr);
  
}











