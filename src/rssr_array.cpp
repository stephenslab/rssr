#include <RcppEigen.h>
#include "rssr.h"
#include <cstdio>
#include <cmath>
#include <RcppParallel.h>


//' Run RSS with the variational bayes algorithm accelerated with SQUAREM, only returning the lower bound
//' @template rssr
//' @param talpha0 a length p vector specifying the initial value of alpha
//' @param tmu0 a length p vector specifying the initial value of mu
//' @param SiRiSr0 a length p vector specifying the initial value of SiRiSr
Eigen::ArrayXd rss_varbvsr_squarem_iter_array(const c_Matrix_internal SiRiS,
                                              const c_arrayxd_internal sigma_beta,
                                              const c_arrayxd_internal logodds,
                                              const c_arrayxd_internal  betahat,
                                              const c_arrayxd_internal se,
                                              const c_arrayxd_internal talpha0,
                                              const c_arrayxd_internal tmu0,
                                              const c_arrayxd_internal tSiRiSr0,
                                              double tolerance,
                                              int itermax,
                                              Rcpp::LogicalVector lnz_tol){
  
  
  //This function implements RSS with variational bayes and the SQUAREM algorithm.
  
  bool lnztol=lnz_tol[0];
  const size_t p = betahat.size();
  const size_t gsize=sigma_beta.size();
  if(gsize!=logodds.size()){
    Rcpp::stop("sigma_beta must be same size as logodds");
  }
  
  
  Eigen::ArrayXXd alpha(p,gsize);
  Eigen::ArrayXXd mu(p,gsize);
  Eigen::ArrayXXd SiRiSr(p,gsize);
  for(int i=0;i<gsize;i++){
    alpha.col(i)=talpha0;
    mu.col(i)=tmu0;
    SiRiSr.col(i)=tSiRiSr0;
  }

  Eigen::ArrayXd lnZ(gsize);
  
  Eigen::ArrayXXd alpha0=alpha;
  Eigen::ArrayXXd mu0=mu;
  Eigen::ArrayXXd SiRiSr0=SiRiSr;
  
  Eigen::ArrayXXd alpha1=alpha;
  Eigen::ArrayXXd mu1=mu;
  Eigen::ArrayXXd SiRiSr1=SiRiSr;

  
  Eigen::ArrayXXd alpha_r(p,gsize);
  Eigen::ArrayXXd alpha_v(p,gsize);
  
  Eigen::ArrayXXd mu_v(p,gsize);
  Eigen::ArrayXXd mu_r(p,gsize);
  
  Eigen::ArrayXd sesquare =se*se;
  Eigen::ArrayXd  q= betahat/sesquare;
  Eigen::ArrayXXd  s(p,gsize);
  for(int i=0; i<gsize;i++){
    s.col(i)= (sesquare*(sigma_beta.coeff(i)*sigma_beta.coeff(i)))/(sesquare+(sigma_beta.coeff(i)*sigma_beta.coeff(i)));
    lnZ[i]=calculate_lnZ(q,alpha.col(i)*mu.col(i),SiRiSr.col(i),logodds.coeff(i),sesquare,alpha.col(i),mu.col(i),s.col(i),sigma_beta.coeff(i));
  }
  


  Eigen::ArrayXd mtp(gsize);  
  size_t iter=0;
  double max_err=1;
  Eigen::ArrayXd lnZ0=lnZ;

//  double rel_l0=0;
  Eigen::ArrayXi rel_li(gsize);
  rel_li.setZero();
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
    for(int i=0; i<gsize;i++){
      mtp.coeffRef(i)= -sqrt((alpha_r.square()).sum()+(mu_r.square()).sum())/sqrt((alpha_v.square()).sum()+(mu_v.square()).sum());
    }
    for(int i=0; i<gsize; i++){
      if(mtp.coeff(i) >=-1){
        
      }else{
        alpha.col(i)=alpha0.col(i)-2*mtp.coeff(i)*alpha_r.col(i)+(mtp.coeff(i)*mtp.coeff(i))*alpha_v.col(i);
        mu.col(i)=mu0.col(i)-2*mtp.coeff(i)*mu_r.col(i)+(mtp.coeff(i)*mtp.coeff(i))*mu_v.col(i);
        SiRiSr.col(i)=SiRiS*(alpha.col(i)*mu.col(i)).matrix();
      }
    }
    
    rss_varbvsr_iter(SiRiS,sigma_beta,logodds,betahat,se,alpha,mu,SiRiSr,reverse);
    for(int i=0;i<gsize;i++){
      lnZ.coeffRef(i)=  calculate_lnZ(q,alpha.col(i)*mu.col(i),SiRiSr.col(i),logodds.coeff(i),sesquare,alpha.col(i),mu.col(i),s.col(i),sigma_beta.coeff(i));
    }
    // if(!std::isfinite(lnZ)){
    //   Rcpp::stop("lnZ isn't finite!");
    // }
    for(int i=0; i<gsize; i++){
      double tmtp = mtp.coeff(i);
      double tlnZ = lnZ.coeff(i);
      double tlnZ0 = lnZ0.coeff(i);
      if((tmtp<(-1)) && (tlnZ < tlnZ0)){
        size_t num_bt=0;
        while((tlnZ<tlnZ0) &&(num_bt < 10)){
          tmtp = 0.5*(tmtp-1);
          alpha.col(i) = alpha0.col(i)-2*tmtp*alpha_r.col(i)+(tmtp*tmtp)*alpha_v.col(i);
          mu.col(i) = mu0.col(i)-2*tmtp*mu_r.col(i)+(tmtp*tmtp)*mu_v.col(i);
          SiRiSr.col(i) = SiRiS*(alpha.col(i)*mu.col(i)).matrix();
          rss_varbvsr_iter(SiRiS,sigma_beta.coeff(i),logodds.coeff(i),betahat,se,alpha.col(i),mu.col(i),SiRiSr.col(i),reverse);
          lnZ.coeffRef(i)=  calculate_lnZ(q,alpha.col(i)*mu.col(i),SiRiSr.col(i),logodds.coeff(i),sesquare,alpha.col(i),mu.col(i),s.col(i),sigma_beta.coeff(i));
          tlnZ=lnZ.coeff(i);
          num_bt=num_bt+1;
          //          lnZ=calculate_lnZ(q,alpha.col(i)*mu.col(i),SiRiSr.col(i),logodds.coeff(i),sesquare,alpha.col(i),mu.col(i),s.col(i),sigma_beta.coeff(i));          
        }
      }
      rel_li.coeffRef(i)=rel_err(lnZ.coeff(i),lnZ0.coeff(i));
    }

    
    if(lnztol){
      max_err=rel_li.maxCoeff();
    }else{
      max_err=find_maxerr(alpha,alpha0,alpha*mu,alpha0*mu0);
    }
    
    if(iter>itermax){
      printf("Maximum iteration number reached: %+0.2d \n",(int)iter);
      printf("The log variational lower bound of the last step increased by %+0.2e\n",(lnZ-lnZ0).maxCoeff());
      break;
    }
    iter=iter+1;
  }
  return lnZ;
}
