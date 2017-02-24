#include <RcppEigen.h>
#include "rssvarbvsr.hpp"
#include "kl.hpp"
#include <cstdio>
#include "tbb/tbb.h"

using namespace tbb;

//' Run RSS with the variational bayes algorithm accelerated with SQUAREM, only returning the lower bound
//' @template rssr
//' @param talpha0 a length p vector specifying the initial value of alpha
//' @param tmu0 a length p vector specifying the initial value of mu
//' @param SiRiSr0 a length p vector specifying the initial value of SiRiSr
//' @useDynLib rssr
//[[Rcpp::export]]
double rss_varbvsr_squarem_grid(const Eigen::MappedSparseMatrix<double> SiRiS,
                                const double sigma_beta,
                                const double logodds,
                                const Eigen::Map<Eigen::ArrayXd> betahat,
                                const Eigen::Map<Eigen::ArrayXd> se,
                                const Eigen::Map<Eigen::ArrayXd> talpha0,
                                const Eigen::Map<Eigen::ArrayXd> tmu0,
                                const Eigen::Map<Eigen::ArrayXd> tSiRiSr0,
                                double tolerance,
                                int itermax,
                                Rcpp::LogicalVector lnz_tol){
  
  
  //This function implements RSS with variational bayes and the SQUAREM algorithm.
  
  bool lnztol=lnz_tol[0];
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
  
  
  lnZ=calculate_lnZ(q,alpha*mu,SiRiSr,logodds,sesquare,alpha,mu,s,sigma_beta);
  double mtp;  
  size_t iter=0;
  double max_err=1;
  double lnZ0=lnZ;
  double lnZ00=lnZ0;
  double rel_l0=0;
  double rel_li=0;
  while(max_err>tolerance){
    lnZ0=lnZ;
    alpha0=alpha;
    mu0=mu; 
    bool reverse = iter%2!=0;
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
    rel_l0=rel_err(lnZ00,lnZ);
    rel_li=rel_err(lnZ,lnZ0);
    if(lnztol){
      max_err=rel_li;
    }else{
      max_err=find_maxerr(alpha,alpha0,alpha*mu,alpha0*mu0);
    }

    if(iter>itermax){
      printf("Maximum iteration number reached: %+0.2d \n",(int)iter);
      printf("The log variational lower bound of the last step increased by %+0.2e\n",lnZ-lnZ0);
      break;
    }
    iter=iter+1;
  }
  return lnZ;
}

//[[Rcpp::export]]
Rcpp::NumericMatrix grid_search_rss_varbvsr(
    const Eigen::MappedSparseMatrix<double> SiRiS,
    const Eigen::Map<Eigen::ArrayXd> sigma_beta,
    const Eigen::Map<Eigen::ArrayXd> logodds,
    const Eigen::Map<Eigen::ArrayXd> betahat,
    const Eigen::Map<Eigen::ArrayXd> se,
    const Eigen::Map<Eigen::ArrayXd> talpha0,
    const Eigen::Map<Eigen::ArrayXd> tmu0,
    const Eigen::Map<Eigen::ArrayXd> tSiRiSr0,
    double tolerance,
    int itermax,
    Rcpp::LogicalVector verbose,
    Rcpp::LogicalVector lnz_tol){
  
  
  size_t sigb_size= sigma_beta.size();
  size_t logodds_size=logodds.size();
  size_t tot_size=sigb_size*logodds_size;

  Rcpp::NumericMatrix nlzmat(logodds_size,sigb_size);
  parallel_for(blocked_range<size_t>(0,tot_size),
               [&](const blocked_range<size_t>& r){
                 for(size_t t=r.begin(); t!=r.end(); t++){
                   size_t i=t/logodds_size;
                   size_t j=t%logodds_size;
                   nlzmat(i,j)=rss_varbvsr_squarem_grid(SiRiS,
                          sigma_beta(j),
                          logodds(i),
                          betahat,
                          se,
                          talpha0,
                          tmu0,
                          tSiRiSr0,
                          tolerance,
                          itermax,
                          lnz_tol);}});
  return(nlzmat);
}



