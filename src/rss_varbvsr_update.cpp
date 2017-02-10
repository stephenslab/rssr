// For a description of this C code, see rss_varbvsr_update.m.
#include <RcppEigen.h>
#include "rssvarbvsr.hpp"
#include "kl.hpp"
#include <cstdio>
                     
void rss_varbvsr_iter(const Eigen::SparseMatrix<double> SiRiS,
                      const double sigma_beta,
                      const double logodds,
                      const Eigen::ArrayXd betahat,
                      const Eigen::ArrayXd se,
                      Eigen::ArrayXd &alpha,
                      Eigen::ArrayXd &mu,
                      Eigen::ArrayXd &SiRiSr,
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
    double SiRiSr_snp = SiRiSr(i);
    
    // Perform the mean-field variational update.
    rss_varbvsr_update(betahat(i), se(i), sigma_beta, 
                       SiRiS_snp, SiRiSr, SiRiSr_snp, 
                       logodds, alpha(i), mu(i));
  }
}



//' Run RSS with the variational bayes algorithm accelerated with SQUAREM
//' @template rssr
//' @param talpha0 a length p vector specifying the initial value of alpha
//' @param tmu0 a length p vector specifying the initial value of mu
//' @param SiRiSr0 a length p vector specifying the initial value of SiRiSr
//' @useDynLib rssr
//[[Rcpp::export]]
Rcpp::List rss_varbvsr_squarem(const Eigen::SparseMatrix<double> &SiRiS,
                                   const double sigma_beta,
                                   const double logodds,
                                   const Eigen::ArrayXd betahat,
                                   const Eigen::ArrayXd se,
                                   const Eigen::ArrayXd &talpha0,
                                   const Eigen::ArrayXd &tmu0,
                                   const Eigen::ArrayXd &tSiRiSr0,
                                   double tolerance){
  
  
  //This function implements RSS with variational bayes and the SQUAREM algorithm.
  
  
  const size_t p = betahat.size();
  
  Eigen::ArrayXd alpha=talpha0;
  Eigen::ArrayXd mu=tmu0;
  Eigen::ArrayXd SiRiSr=tSiRiSr0;
  double lnZ=log(0);
  Eigen::MatrixXd params0(p,2);
  Eigen::MatrixXd params(p,2);
  
  Eigen::ArrayXd alpha0=alpha;
  Eigen::ArrayXd mu0=mu;
  Eigen::ArrayXd SiRiSr0=SiRiSr;
  
  Eigen::ArrayXd alpha1=alpha;
  Eigen::ArrayXd mu1=mu;
  Eigen::ArrayXd SiRiSr1=SiRiSr;
  
  Eigen::ArrayXd alpha3=alpha;
  Eigen::ArrayXd mu3=mu;
  Eigen::ArrayXd SiRiSr3=SiRiSr;
  
  Eigen::ArrayXd alpha_r(p);
  Eigen::ArrayXd alpha_v(p);
  
  Eigen::ArrayXd mu_v(p);
  Eigen::ArrayXd mu_r(p);

  Eigen::ArrayXd sesquare =se*se;
  Eigen::ArrayXd  q= betahat/sesquare;
  Eigen::ArrayXd  s= (sesquare*(sigma_beta*sigma_beta))/(sesquare+(sigma_beta*sigma_beta));
  
  
  
  double mtp;  
  size_t iter=0;
  double max_err=1;
  double lnZ0=0;
  bool  verbose=true;
  while(max_err>tolerance){
    lnZ0=lnZ;
    alpha0=alpha;
    mu0=mu;
    bool reverse = iter%2==0;
    rss_varbvsr_iter(SiRiS,sigma_beta,logodds,betahat,se,alpha,mu,SiRiSr,reverse);
    alpha1=alpha;
    mu1=mu;
    SiRiSr1=SiRiSr;
    
    rss_varbvsr_iter(SiRiS,sigma_beta,logodds,betahat,se,alpha,mu,SiRiSr,reverse);
    
    alpha_r=alpha1-alpha0;
    alpha_v=(alpha-alpha1)-alpha_r;
      
    mu_r   = mu1-mu0;
    mu_v   = mu-mu1-mu_r;
    
    mtp= -sqrt((alpha_r.square()).sum()+(mu_r.square()).sum())/sqrt((alpha_v.square()).sum()+(mu_v.square()).sum());

    if(mtp >=-1){
      
    }else{
      alpha=alpha0-2*mtp*alpha_r+(mtp*mtp)*alpha_v;
      mu=mu0-2*mtp*mu_r+(mtp*mtp)*mu_v;
      SiRiSr=SiRiS*(alpha*mu).matrix();
    }
    
    rss_varbvsr_iter(SiRiS,sigma_beta,logodds,betahat,se,alpha,mu,SiRiSr,reverse);
    lnZ=  calculate_lnZ(q,alpha*mu,SiRiSr,logodds,sesquare,alpha,mu,s,sigma_beta);
    if((mtp<(-1)) && (lnZ < lnZ0)){
      size_t num_bt=0;
      while((lnZ<lnZ0) &&(num_bt < 10)){
        mtp = 0.5*(mtp-1);
        alpha = alpha0-2*mtp*alpha_r+(mtp*mtp)*alpha_v;
        mu = mu0-2*mtp*mu_r+(mtp*mtp)*mu_v;
        SiRiSr = SiRiS*(alpha*mu).matrix();
        rss_varbvsr_iter(SiRiS,sigma_beta,logodds,betahat,se,alpha,mu,SiRiSr,reverse);
        lnZ=calculate_lnZ(q,alpha*mu,SiRiSr,logodds,sesquare,alpha,mu,s,sigma_beta);          
        num_bt=num_bt+1;
      }
    }
    
    max_err=find_maxerr(alpha,alpha0,alpha*mu,alpha0*mu0);
    if(verbose){
      double absr=(alpha*mu).abs().maxCoeff();
      int asum=round(alpha.sum());
      printf("%4d %+13.6e %0.1e %4d %0.2f %5.2f\n",(int)iter,lnZ,max_err,(int) asum,absr,sigma_beta*sigma_beta);
    }
    iter=iter+1;
  }
  return Rcpp::List::create(Rcpp::Named("alpha")=alpha,
                     Rcpp::Named("mu")=mu,
                     Rcpp::Named("SiRiSr")=SiRiSr,
                     Rcpp::Named("max_err")=max_err,
                     Rcpp::Named("lnZ")=lnZ,
                     Rcpp::Named("iter")=iter);
}

//[[Rcpp::export]]
Eigen::MatrixXd rss_varbvsr_naive (const Eigen::SparseMatrix<double> &SiRiS,
                                  const double sigma_beta,
                                  const double logodds,
                                  const Eigen::ArrayXd betahat,
                                  const Eigen::ArrayXd se,
                                  const Eigen::ArrayXd &alpha0,
                                  const Eigen::ArrayXd &mu0,
                                  const Eigen::ArrayXd &SiRiSr0,
                                  double tolerance){
  
  // This is the variational implementation in pure C++, it does not have the SQUAREM update.
  
  
  // Get the number of SNPs (p) and coordinate ascent updates (m).
  const size_t p = betahat.size();
  Eigen::ArrayXd alpha=alpha0;
  Eigen::ArrayXd mu=mu0;
  Eigen::ArrayXd SiRiSr=SiRiSr0;
  double lnZ=log(0);
  Eigen::MatrixXd params0(p,2);
  Eigen::MatrixXd params(p,2);
  
  Eigen::ArrayXd talpha=alpha;
  Eigen::ArrayXd tmu=mu;
  
  
  size_t iter=0;
  //  Initialize outputs.
  double max_err=1;
  while(max_err>tolerance){
    talpha=alpha;
    tmu=mu;
    bool reverse = iter%2==0;
    rss_varbvsr_iter(SiRiS,sigma_beta,logodds,betahat,se,alpha,mu,SiRiSr,reverse);
    max_err=find_maxerr(alpha,talpha,alpha*mu,talpha*tmu);
  }
  Eigen::MatrixXd retmat(p,3);
  retmat << alpha,mu,SiRiSr;
  return(retmat);
}
