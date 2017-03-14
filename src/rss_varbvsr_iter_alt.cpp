#include <RcppEigen.h>
#include "kl.hpp"
#include "rss_varbvsr_iter_alt.hpp"
#include "sigmoid.hpp"
#include "sparse_matrix.hpp"
#include <math.h> 

//[[Rcpp::interfaces(r,cpp)]]




//
//
//[[Rcpp::export]]
Rcpp::List rss_varbvsr_naive_alt (const Eigen::MappedSparseMatrix<double> &SiRiS,
                                  const double sigma_beta,
                                  const double logodds,
                                  const Eigen::Map<Eigen::ArrayXd> betahat,
                                  const Eigen::Map<Eigen::ArrayXd> se,
                                  const Eigen::ArrayXd &talpha0,
                                  const Eigen::ArrayXd &tmu0,
                                  const Eigen::ArrayXd &tSiRiSr0,
                                  double tolerance,
                                  int itermax,
                                  Rcpp::LogicalVector verbose,
                                  Rcpp::LogicalVector lnz_tol){
  
  // This is the variational implementation in pure C++, it does not have the SQUAREM update.
  
  
  
  
  // Get the number of SNPs (p) and coordinate ascent updates (m).
  const size_t p = betahat.size();
  
  
  Eigen::ArrayXd alpha=talpha0;
  Eigen::ArrayXd mu=tmu0;
  Eigen::ArrayXd SiRiSr=tSiRiSr0;
  double lnZ=0;
  
  
  Eigen::ArrayXd talpha=alpha;
  Eigen::ArrayXd tmu=mu;
  Eigen::ArrayXd sesquare=se*se;
  
  Eigen::ArrayXd  q= betahat/sesquare;
  Eigen::ArrayXd  s= (sesquare*(sigma_beta*sigma_beta))/(sesquare+(sigma_beta*sigma_beta));
  bool verbosev=verbose[0];
  bool lnztol=lnz_tol[0];
  size_t iter=0;
  //  Initialize outputs.
  double max_err=1;
  lnZ=calculate_lnZ(q,alpha*mu,SiRiSr.array(),logodds,sesquare,alpha,mu,s,sigma_beta);
  
  double lnZ0=lnZ;
  double lnZ00=lnZ0;
  double rel_l0=0;
  double rel_li=0;
  while(max_err>tolerance){
    talpha=alpha;
    tmu=mu;
    bool reverse = iter%2!=0;
    rss_varbvsr_iter_alt(  SiRiS,
                           sigma_beta,
                           logodds,
                           betahat,
                           se,
                           alpha,
                           mu,
                           SiRiSr,
                           reverse);
    rel_li=rel_err(lnZ,lnZ0);
    if(lnztol){
      max_err=rel_li;
    }else{
      max_err=find_maxerr(alpha,talpha,alpha*mu,talpha*tmu);
    }
    
    iter=iter+1;
  }
  lnZ=  calculate_lnZ(q,alpha*mu,SiRiSr.array(),logodds,sesquare,alpha,mu,s,sigma_beta);
  
  return Rcpp::List::create(Rcpp::Named("alpha")=alpha,
                            Rcpp::Named("mu")=mu,
                            Rcpp::Named("SiRiSr")=SiRiSr,
                            Rcpp::Named("max_err")=max_err,
                            Rcpp::Named("lnZ")=lnZ,
                            Rcpp::Named("iter")=iter);
  
}


void compute_mu(const double betahat,
                const double se_square,
                const double sigma_square,
                const double alpha,
                double &mu,
                const double SiRiSr_snp){
  
  
  
  // Compute the variational estimate of the posterior variance.
  
  mu= sigma_square * (betahat / se_square + (alpha*mu) / se_square - SiRiSr_snp);
}

void compute_alpha(const double sigma_square,
                   const double sigma_beta,
                   const double logodds,
                   const double mu,
                   double &alpha) {
  
  double sigma_beta_square = sigma_beta * sigma_beta;
  
  // Compute the variational estimate of the posterior variance.

  
  // Update the variational estimate of the posterior inclusion probability.
  alpha=sigmoid(logodds + 0.5 * (log(sigma_square/(sigma_beta_square)) + (mu*mu)/sigma_square));
}

// void compute_SiRiSr(const Eigen::SparseVector<double> &SiRiS_snp,
//                     const  double r0,
//                     const double r_new,
//                     Eigen::ArrayXd &SiRiSr) {
//   SiRiSr+=(SiRiS_snp*(r_new-r0));
// 
// }




void rss_varbvsr_iter_alt(const Eigen::SparseMatrix<double> &SiRiS,
                          const double sigma_beta,
                          const double logodds,
                          const Eigen::Map<Eigen::ArrayXd> betahat,
                          const Eigen::Map<Eigen::ArrayXd> se,
                          Eigen::ArrayXd &alpha,
                          Eigen::ArrayXd &mu,
                          Eigen::ArrayXd &SiRiSr,
                          bool reverse){
  
  
  // Get the number of SNPs (p) and coordinate ascent updates (m).
  const size_t p = betahat.size();
  
  // Initialize outputs.
  
  // Store a single column of matrix inv(S)*R*inv(S).
  // Eigen::ArrayXd alpha0=alpha;
  // Eigen::ArrayXd mu0=mu;
  //  Eigen::ArrayXd SiRiSr0=SiRiSr;
  double sigma_beta_square = sigma_beta * sigma_beta;
  Eigen::ArrayXd sesquare=se*se;
  Eigen::ArrayXd sigma_square = (sesquare * sigma_beta_square) / (sesquare + sigma_beta_square);
  
  size_t i=0;
  
  double r0=0;
  double r_new=0;
  //  double sigma_beta_square = sigma_beta * sigma_beta;
  double *talpha=NULL;
  double *tmu=NULL;
  const double *tbetahat=NULL;
  const double *tsesquare=NULL;
  double *tSiRiSr_snp=NULL;
  const double *tsigma_square=NULL;
  Eigen::ArrayXd  SiRiS_snp(p);
  Eigen::VectorXd  SiRiS_snp_v(p);
  if(!reverse){
    for(size_t i=0; i<p; i++){
      
      talpha=&alpha.coeffRef(i);
      tmu = &mu.coeffRef(i);
      tsigma_square=&sigma_square.coeffRef(i);
      r0=(*talpha)*(*tmu);
      tbetahat= &betahat.coeffRef(i);
      tsesquare=&sesquare.coeffRef(i);
      tSiRiSr_snp=&SiRiSr.coeffRef(i);
      
      
      *tmu= *tsigma_square * ((*tbetahat) / (*tsesquare) + ((r0) / (*tsesquare) - *tSiRiSr_snp));
      //      compute_mu(betahat(i),sesquare(i),sigma_square(i),alpha(i),mu(i),SiRiSr(i));
      (*talpha)= sigmoid(logodds + 0.5 * (log(*tsigma_square/sigma_beta_square) + ((*tmu * *tmu)/(*tsigma_square))));
      r_new=*talpha* (*tmu);
      SiRiS_snp_v=SiRiS.col(i);
      SiRiSr+=(SiRiS_snp_v.array()*(r_new-r0));
      // compute_SiRiSr(SiRiS.col(i),r0,r_new,SiRiSr);
    }
  } else{
    for(size_t i=p-1; i>0; i--){
      
      talpha=&alpha.coeffRef(i);
      tmu = &mu.coeffRef(i);
      tsigma_square=&sigma_square.coeffRef(i);
      r0=(*talpha)*(*tmu);
      tbetahat= &betahat.coeffRef(i);
      tsesquare=&sesquare.coeffRef(i);
      tSiRiSr_snp=&SiRiSr.coeffRef(i);
      
      
      *tmu= *tsigma_square * ((*tbetahat) / (*tsesquare) + ((r0) / (*tsesquare) - *tSiRiSr_snp));
      //      compute_mu(betahat(i),sesquare(i),sigma_square(i),alpha(i),mu(i),SiRiSr(i));
      (*talpha)= sigmoid(logodds + 0.5 * (log(*tsigma_square/sigma_beta_square) + ((*tmu * *tmu)/(*tsigma_square))));
      r_new=*talpha* (*tmu);
      SiRiS_snp_v=SiRiS.col(i);
      SiRiSr+=(SiRiS_snp_v.array()*(r_new-r0));
    }
  }
  
  
}








//[[Rcpp::export]]
Rcpp::List wrap_rss_varbvsr_iter_alt(const Eigen::SparseMatrix<double> &SiRiS,
                                     const double sigma_beta,
                                     const double logodds,
                                     const Eigen::Map<Eigen::ArrayXd> betahat,
                                     const Eigen::Map<Eigen::ArrayXd> se,
                                     Eigen::ArrayXd &alpha,
                                     Eigen::ArrayXd &mu,
                                     Eigen::ArrayXd &SiRiSr,
                                     bool reverse){

//  Eigen::ArrayXd sesquare=se*se;
rss_varbvsr_iter_alt(  SiRiS,
                       sigma_beta,
                       logodds,
                       betahat,
                       se,
                       alpha,
                       mu,
                       SiRiSr,
                       reverse);
  return Rcpp::List::create(Rcpp::Named("alpha1")=alpha,
                            Rcpp::Named("mu1")=mu,
                            Rcpp::Named("SiRiSr")=SiRiSr);
}

// [[Rcpp::export]]
double wrap_compute_alpha(const double sigma_square,
                          const double sigma_beta,
                          const double logodds,
                          const double mu,
                          const double alpha){
  double talpha=alpha;
  compute_alpha(sigma_square,sigma_beta,logodds,mu,talpha);
  return(talpha);
}

//[[Rcpp::export]]
double wrap_compute_mu(const double betahat,
                       const double se_square,
                       const double sigma_square,
                       const double alpha,
                       const double mu,
                       const double SiRiSr_snp){
  double tmu=mu;
  compute_mu(betahat,se_square,sigma_square,alpha,tmu,SiRiSr_snp);
  return(tmu);
}
// 

// Eigen::ArrayXd wrap_compute_SiRiSr(const Eigen::SparseMatrix<double> SiRiS_snp,
//                                    const  double r0,
//                                    const double r_new,
//                                    const Eigen::ArrayXd SiRiSr){
//   Eigen::VectorXd tSiRiSr=SiRiSr;
//   compute_SiRiSr(SiRiS_snp.col(0),r0,r_new,tSiRiSr);
//   return(tSiRiSr.array());
// }

  

    

