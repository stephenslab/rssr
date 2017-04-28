#include <RcppEigen.h>
#include "rssr.h"
#include <cstdio>
#include <cmath>
#include <RcppParallel.h>
//[[Rcpp::depends(RcppParallel)]]






//' Run RSS with the variational bayes algorithm accelerated with SQUAREM, only returning the lower bound
//' @template rssr
//' @param talpha0 a length p vector specifying the initial value of alpha
//' @param tmu0 a length p vector specifying the initial value of mu
//' @param SiRiSr0 a length p vector specifying the initial value of SiRiSr
double rss_varbvsr_squarem_iter(const c_sparseMatrix_internal SiRiS,
                                   const double sigma_beta,
                                   const double logodds,
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
    // if(!std::isfinite(lnZ)){
    //   Rcpp::stop("lnZ isn't finite!");
    // }
    
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








double rss_varbvsr_naive_iter (const Matrix_external SiRiS,
                              const double sigma_beta,
                              const double logodds,
                              const arrayxd_external betahat,
                              const arrayxd_external se,
                              const arrayxd_external talpha0,
                              const arrayxd_external tmu0,
                              const arrayxd_external tSiRiSr0,
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
  lnZ=calculate_lnZ(q,alpha*mu,SiRiSr,logodds,sesquare,alpha,mu,s,sigma_beta);
  
  double lnZ0=lnZ;
  double lnZ00=lnZ0;
  double rel_l0=0;
  double rel_li=0;
  while(max_err>tolerance){
    talpha=alpha;
    tmu=mu;
    lnZ0=lnZ;
    bool reverse = iter%2!=0;
    rss_varbvsr_iter(SiRiS,sigma_beta,logodds,betahat,se,alpha,mu,SiRiSr,reverse);
    lnZ=  calculate_lnZ(q,alpha*mu,SiRiSr,logodds,sesquare,alpha,mu,s,sigma_beta);
    rel_li=rel_err(lnZ,lnZ0);
    if(lnztol){
      max_err=rel_li;
    }else{
      max_err=find_maxerr(alpha,talpha,alpha*mu,talpha*tmu);
    }
    
    iter=iter+1;
  }
  lnZ=  calculate_lnZ(q,alpha*mu,SiRiSr,logodds,sesquare,alpha,mu,s,sigma_beta);
  
  return lnZ;
  
}





//' Run RSS with the variational bayes algorithm accelerated with SQUAREM, only returning the lower bound
//' @template rssr
//' @param talpha0 a length p vector specifying the initial value of alpha
//' @param tmu0 a length p vector specifying the initial value of mu
//' @param SiRiSr0 a length p vector specifying the initial value of SiRiSr
double rss_varbvsr_squarem_iter(const c_Matrix_internal SiRiS,
                                const double sigma_beta,
                                const double logodds,
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
  

  Eigen::VectorXd r=alpha*mu;
  
  Eigen::ArrayXd alpha_r(p);
  Eigen::ArrayXd alpha_v(p);
  
  Eigen::ArrayXd mu_v(p);
  Eigen::ArrayXd mu_r(p);
  double sigma_beta_square=sigma_beta*sigma_beta;
  Eigen::ArrayXd sesquare =se*se;
  Eigen::ArrayXd  q= betahat/sesquare;
  Eigen::ArrayXd  s= (sesquare*(sigma_beta*sigma_beta))/(sesquare+(sigma_beta*sigma_beta));
  Eigen::ArrayXd ssrat((s/sigma_beta_square).log());
  if(r.hasNaN()){
    Rcpp::Rcerr<<"In iteration iter 0(0)"<<std::endl;
    Rcpp::stop("alpha*mu is not finite!");
  }

  lnZ=calculate_lnZ(q,r,SiRiSr,logodds,sesquare,alpha,mu,s,sigma_beta);
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
    rss_varbvsr_iter(SiRiS,sigma_beta_square,s,logodds,betahat,sesquare,ssrat,alpha,mu,SiRiSr,reverse);
//    r=alpha*mu;
    alpha1=alpha;
    mu1=mu;
    SiRiSr1=SiRiSr;
    rss_varbvsr_iter(SiRiS,sigma_beta_square,s,logodds,betahat,sesquare,ssrat,alpha,mu,SiRiSr,reverse);
    
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
    
    rss_varbvsr_iter(SiRiS,sigma_beta_square,s,logodds,betahat,sesquare,ssrat,alpha,mu,SiRiSr,reverse);
    r=alpha*mu;

    // if(r.hasNaN()){
    //   Rcpp::Rcerr<<"In iteration iter: "<<iter<<std::endl;
    //   Rcpp::Rcerr<<"alpha.min()"<<alpha.minCoeff()<<std::endl;
    //   Rcpp::Rcerr<<"alpha.max()"<<alpha.maxCoeff()<<std::endl;
    //   Rcpp::Rcerr<<"mu.min()"<<mu.minCoeff()<<std::endl;
    //   Rcpp::Rcerr<<"mu.min()"<<mu.maxCoeff()<<std::endl;
    //   Rcpp::stop("alpha*mu is not finite!");
    // }
  
  
    lnZ=  calculate_lnZ(q,r,SiRiSr,logodds,sesquare,alpha,mu,s,sigma_beta);
    // if(!std::isfinite(lnZ)){
    //   Rcpp::stop("lnZ isn't finite!");
    // }
    
    if((mtp<(-1)) && (lnZ < lnZ0)){
      size_t num_bt=0;
      while((lnZ<lnZ0) &&(num_bt < 10)){
        mtp = 0.5*(mtp-1);
        alpha = alpha0-2*mtp*alpha_r+(mtp*mtp)*alpha_v;
        mu = mu0-2*mtp*mu_r+(mtp*mtp)*mu_v;
        SiRiSr = SiRiS*(alpha*mu).matrix();
        rss_varbvsr_iter(SiRiS,sigma_beta_square,s,logodds,betahat,se,ssrat,alpha,mu,SiRiSr,reverse);
        r=alpha*mu;
        // if(r.hasNaN()){
        //   Rcpp::Rcerr<<"In iteration iter: "<<iter<<" (num_bt: "<<num_bt<<")"<<std::endl;
        //   Rcpp::stop("alpha*mu is not finite!");
        // }
        lnZ=calculate_lnZ(q,r,SiRiSr,logodds,sesquare,alpha,mu,s,sigma_beta);          
        num_bt=num_bt+1;
      }
    }
    rel_l0=rel_err(lnZ00,lnZ);
    rel_li=rel_err(lnZ,lnZ0);
    
    if(lnztol){
      max_err=rel_li;
    }else{
      max_err=find_maxerr(alpha,alpha0,r,alpha0*mu0);
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





//' Run RSS with the variational bayes algorithm accelerated with SQUAREM, only returning the lower bound
//' @template rssr
//' @param talpha0 a length p vector specifying the initial value of alpha
//' @param tmu0 a length p vector specifying the initial value of mu
//' @param SiRiSr0 a length p vector specifying the initial value of SiRiSr
Eigen::MatrixXd rss_varbvsr_squarem_iter_trace(const c_Matrix_internal SiRiS,
                                const double sigma_beta,
                                const double logodds,
                                const c_arrayxd_internal  betahat,
                                const c_arrayxd_internal se,
                                const c_arrayxd_internal talpha0,
                                const c_arrayxd_internal tmu0,
                                const c_arrayxd_internal tSiRiSr0,
                                double tolerance,
                                int itermax,
                                Rcpp::LogicalVector lnz_tol){
  

  Eigen::MatrixXd retmat(itermax,3);
  retmat.setZero();
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
    if(!std::isfinite(lnZ)){
      Rcpp::stop("lnZ isn't finite!");
    }
    
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
    if(alpha.mean()>=1){
      Rcpp::stop("Mean of alpha can't be >= 1");
    }
    double tsigb_square = (alpha.matrix().dot((s+mu.square()).matrix()))/(alpha).sum();
    retmat.coeffRef(iter,0)=lnZ;
    retmat.coeffRef(iter,1)=alpha.mean();
    retmat.coeffRef(iter,2)=sqrt(tsigb_square);
    iter=iter+1;
//    retZ.push_back(lnZ);
  }
  return retmat.topRows(iter);
}








double rss_varbvsr_squarem_iter_fit_logodds(const c_sparseMatrix_internal SiRiS,
                                const double sigma_beta,
                                double &logodds,
                                const c_arrayxd_internal betahat,
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
    logodds=update_logodds(alpha);
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





// 
// 
// Rcpp::List wrap_rss_varbvsr_iter_grid(const sparseMatrix_external SiRiS,
//                                       const arrayxd_external sigma_beta,
//                                       const arrayxd_external logodds,
//                                       const arrayxd_external betahat,
//                                       const arrayxd_external se,
//                                       const arrayxd_external alpha,
//                                       const arrayxd_external mu,
//                                       const arrayxd_external SiRiSr,
//                                       Rcpp::LogicalVector reverse){
//   
//   int p = betahat.size();
//   int tot_size=sigma_beta.size();
//   Eigen::ArrayXd sesquare=se*se;
//   Eigen::ArrayXXd talpha(p,tot_size);
//   talpha.colwise()=alpha;
//   
//   Eigen::ArrayXXd tmu(p,tot_size);
//   tmu.colwise()=mu;
//   Eigen::ArrayXXd tSiRiSr(p,tot_size);
//   tSiRiSr.colwise()=SiRiSr;
//   
//   Eigen::ArrayXXd s(p,tot_size);
//   for(int i=0; i<tot_size; i++){
//     s.col(i)=(sesquare*(sigma_beta(i)*sigma_beta(i)))/(sesquare+(sigma_beta(i)*sigma_beta(i)));
//   }
//   
//   rss_varbvsr_iter_alt_grid(SiRiS,
//                             sigma_beta,
//                             logodds,
//                             betahat,
//                             se,
//                             s,
//                             talpha,
//                             tmu,
//                             tSiRiSr,
//                             reverse(0));
//   return Rcpp::List::create(Rcpp::Named("alpha1")=talpha,
//                             Rcpp::Named("mu1")=tmu,
//                             Rcpp::Named("SiRiSr")=tSiRiSr);
// }
// 
// 
// 
// Rcpp::List rss_varbvsr_alt_naive_grid(const c_sparseMatrix_internal SiRiS,
//                                       const arrayxd_external sigma_beta,
//                                       const arrayxd_external logodds,
//                                       const arrayxd_external betahat,
//                                       const arrayxd_external se,
//                                       const arrayxd_external talpha0,
//                                       const arrayxd_external tmu0,
//                                       const arrayxd_external tSiRiSr0,
//                                       double tolerance,
//                                       int itermax,
//                                       bool verbose,
//                                       bool lnztol){
//   
//   
//   //  bool lnztol=lnz_tol[0];
//   using namespace Rcpp;
//   using namespace Eigen;
//   int sigb_size= sigma_beta.size();
//   int logodds_size=logodds.size();
//   int tot_size=sigb_size;
//   if(tot_size!=logodds_size){
//     Rcpp::stop("tot_size!=logodds_size");
//   }
//   int p=betahat.size();
//   
//   ArrayXXd alpha(p,tot_size);
//   alpha.colwise()=talpha0;
//   ArrayXXd mu(p,tot_size);
//   mu.colwise()=tmu0;
//   
//   ArrayXXd SiRiSr(p,tot_size);
//   SiRiSr.colwise()=tSiRiSr0;
//   
//   RowArray lnZ(tot_size);
//   
//   
//   Eigen::ArrayXXd alpha0=alpha;
//   Eigen::ArrayXXd mu0=mu;
//   //  Eigen::ArrayXXd SiRiSr0=SiRiSr;
//   
//   Eigen::ArrayXXd alpha1=alpha;
//   Eigen::ArrayXXd mu1=mu;
//   //  Eigen::ArrayXXd SiRiSr1=SiRiSr;
//   
//   
//   
//   Eigen::ArrayXd sesquare =se*se;
//   Eigen::ArrayXd  q= betahat/sesquare;
//   Eigen::ArrayXXd s(p,tot_size);
//   // Rcout<<"Initializing sigma_square and initial lnZ values"<<std::endl;
//   for(int i=0; i<tot_size; i++){
//     s.col(i)=(sesquare*(sigma_beta(i)*sigma_beta(i)))/(sesquare+(sigma_beta(i)*sigma_beta(i)));
//     lnZ(i)=calculate_lnZ(q,alpha.col(i)*mu.col(i),SiRiSr.col(i),logodds(i),sesquare,alpha.col(i),mu.col(i),s.col(i),sigma_beta(i));
//   }
//   
//   RowArray max_err(tot_size);
//   max_err.setOnes();
//   RowArray lnZ0=lnZ;
//   RowArray rel_li(tot_size);
//   rel_li.setZero();
//   int iter=0;
//   while(max_err.maxCoeff()>tolerance){
//     //    Rcout<<iter<<": worst error is "<<max_err.maxCoeff()<<" tolerance is "<<tolerance<<std::endl;
//     lnZ0=lnZ;
//     alpha0=alpha;
//     mu0=mu;
//     
//     bool reverse = iter%2!=0;
//     
//     rss_varbvsr_iter_alt_grid(SiRiS,sigma_beta,logodds,betahat,se,s,alpha,mu,SiRiSr,reverse);
//     
//     for(int c=0; c<tot_size; c++){
//       lnZ(c)=calculate_lnZ(q,alpha.col(c)*mu.col(c),SiRiSr.col(c),logodds(c),sesquare,alpha.col(c),mu.col(c),s.col(c),sigma_beta(c));
//     }
//     //    rel_l0=rel_err(lnZ00,lnZ);
//     rel_li=(lnZ-lnZ0).abs()/(lnZ.abs()+lnZ0.abs()+double_lim::epsilon());
//     if(lnztol){
//       max_err=rel_li;
//     }else{
//       for(int c=0; c<tot_size; c++){
//         max_err(c)=find_maxerr(alpha.col(c),alpha0.col(c),alpha.col(c)*mu.col(c),alpha0.col(c)*mu0.col(c));
//       }
//     }
//     iter++;
//   }
//   return Rcpp::List::create(Rcpp::Named("alpha")=alpha,
//                             Rcpp::Named("mu")=mu,
//                             Rcpp::Named("SiRiSr")=SiRiSr,
//                             Rcpp::Named("max_err")=max_err,
//                             Rcpp::Named("lnZ")=lnZ,
//                             Rcpp::Named("iter")=iter);
// }








#if RCPP_PARALLEL_USE_TBB

using namespace tbb;


Rcpp::DataFrame grid_rss_varbvsr(
    const c_sparseMatrix_internal SiRiS,
    const c_arrayxd_internal sigma_beta,
    const c_arrayxd_internal logodds,
    const c_arrayxd_internal  betahat,
    const c_arrayxd_internal  se,
    const c_arrayxd_internal talpha0,
    const c_arrayxd_internal tmu0,
    const c_arrayxd_internal tSiRiSr0,
    double tolerance,
    int itermax,
    bool isVerbose,
    bool islnz_tol){
//  std::cout<<"Starting grid_rss_varbvsr (tbb)"<<std::endl;
  
  using namespace Rcpp;
  size_t sigb_size= sigma_beta.size();
  size_t logodds_size=logodds.size();
  size_t tot_size=sigb_size*logodds_size;
  
  Rcpp::NumericVector nlzvec(tot_size);
  Rcpp::NumericVector sigbvec(tot_size);
  Rcpp::NumericVector lovec(tot_size);
  
  Rcpp::LogicalVector verbose(1);
  verbose(0)=isVerbose;
  Rcpp::LogicalVector lnz_tol(1);
  lnz_tol(0)=islnz_tol;
  
  parallel_for(blocked_range<size_t>(0,tot_size),
               [&](const blocked_range<size_t>& r){
                 for(size_t t=r.begin(); t!=r.end(); t++){
                   size_t i=t%logodds_size;
                   size_t j=t/logodds_size;
                   sigbvec(t)=sigma_beta(j);
                   lovec(t)=logodds(i);
                   nlzvec(t)=rss_varbvsr_squarem_iter(SiRiS,
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
  
  //  Rcpp::Rcout<<"mean lnZ is: "<<mean(nlzvec)<<std::endl;
  return(Rcpp::DataFrame::create(_["logodds"]=lovec,
                                 _["sigb"]=sigbvec,
                                 _["lnZ"]=nlzvec));
}




#else


Rcpp::DataFrame grid_rss_varbvsr(
    const c_sparseMatrix_internal SiRiS,
    const c_arrayxd_internal sigma_beta,
    const c_arrayxd_internal logodds,
    const c_arrayxd_internal  betahat,
    const c_arrayxd_internal  se,
    const c_arrayxd_internal talpha0,
    const c_arrayxd_internal tmu0,
    const c_arrayxd_internal tSiRiSr0,
    double tolerance,
    int itermax,
    bool isVerbose,
    bool islnz_tol){
  
  using namespace Rcpp;
  size_t sigb_size= sigma_beta.size();
  size_t logodds_size=logodds.size();
  size_t tot_size=sigb_size*logodds_size;
  
  Rcpp::NumericVector nlzvec(tot_size);
  Rcpp::NumericVector sigbvec(tot_size);
  Rcpp::NumericVector lovec(tot_size);
  
  std::cout<<"Starting grid_rss_varbvsr (serial)"<<std::endl;
  
  Rcpp::LogicalVector verbose(1);
  verbose(0)=isVerbose;
  Rcpp::LogicalVector lnz_tol(1);
  lnz_tol(0)=islnz_tol;
  
  
  for(size_t t=0; t<tot_size; t++){
    size_t i=t%logodds_size;
    size_t j=t/logodds_size;
    sigbvec(t)=sigma_beta(j);
    lovec(t)=logodds(i);
    nlzvec(t)=rss_varbvsr_squarem_iter_sp(SiRiS,
           sigma_beta(j),
           logodds(i),
           betahat,
           se,
           talpha0,
           tmu0,
           tSiRiSr0,
           tolerance,
           itermax,
           lnz_tol);
  }
  
  return(Rcpp::DataFrame::create(_["logodds"]=lovec,
                                 _["sigb"]=sigbvec,
                                 _["lnZ"]=nlzvec));
}

#endif









#if RCPP_PARALLEL_USE_TBB

using namespace tbb;


Rcpp::DataFrame grid_rss_varbvsr(
    const c_Matrix_internal SiRiS,
    const c_arrayxd_internal sigma_beta,
    const c_arrayxd_internal logodds,
    const c_arrayxd_internal  betahat,
    const c_arrayxd_internal  se,
    const c_arrayxd_internal talpha0,
    const c_arrayxd_internal tmu0,
    const c_arrayxd_internal tSiRiSr0,
    double tolerance,
    int itermax,
    bool isVerbose,
    bool islnz_tol){
  
//  std::cout<<"Starting grid_rss_varbvsr (tbb)"<<std::endl;
  
  using namespace Rcpp;
  size_t sigb_size= sigma_beta.size();
  size_t logodds_size=logodds.size();
  size_t tot_size=sigb_size*logodds_size;
  
  Rcpp::NumericVector nlzvec(tot_size);
  Rcpp::NumericVector sigbvec(tot_size);
  Rcpp::NumericVector lovec(tot_size);
  
  Rcpp::LogicalVector verbose(1);
  verbose(0)=isVerbose;
  Rcpp::LogicalVector lnz_tol(1);
  lnz_tol(0)=islnz_tol;
  
  parallel_for(blocked_range<size_t>(0,tot_size),
               [&](const blocked_range<size_t>& r){
                 for(size_t t=r.begin(); t!=r.end(); t++){
                   size_t i=t%logodds_size;
                   size_t j=t/logodds_size;
                   sigbvec(t)=sigma_beta(j);
                   lovec(t)=logodds(i);
                   nlzvec(t)=rss_varbvsr_squarem_iter(SiRiS,
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
  
  //  Rcpp::Rcout<<"mean lnZ is: "<<mean(nlzvec)<<std::endl;
  return(Rcpp::DataFrame::create(_["logodds"]=lovec,
                                 _["sigb"]=sigbvec,
                                 _["lnZ"]=nlzvec));
}




#else


Rcpp::DataFrame grid_rss_varbvsr(
    const c_Matrix_internal SiRiS,
    const c_arrayxd_internal sigma_beta,
    const c_arrayxd_internal logodds,
    const c_arrayxd_internal  betahat,
    const c_arrayxd_internal  se,
    const c_arrayxd_internal talpha0,
    const c_arrayxd_internal tmu0,
    const c_arrayxd_internal tSiRiSr0,
    double tolerance,
    int itermax,
    bool isVerbose,
    bool islnz_tol){
  
  using namespace Rcpp;
  size_t sigb_size= sigma_beta.size();
  size_t logodds_size=logodds.size();
  size_t tot_size=sigb_size*logodds_size;
  
  Rcpp::NumericVector nlzvec(tot_size);
  Rcpp::NumericVector sigbvec(tot_size);
  Rcpp::NumericVector lovec(tot_size);
  
   std::cout<<"Starting grid_rss_varbvsr (serial)"<<std::endl;
  
  Rcpp::LogicalVector verbose(1);
  verbose(0)=isVerbose;
  Rcpp::LogicalVector lnz_tol(1);
  lnz_tol(0)=islnz_tol;
  
  
  for(size_t t=0; t<tot_size; t++){
    size_t i=t%logodds_size;
    size_t j=t/logodds_size;
    sigbvec(t)=sigma_beta(j);
    lovec(t)=logodds(i);
    nlzvec(t)=rss_varbvsr_squarem_iter_sp(SiRiS,
           sigma_beta(j),
           logodds(i),
           betahat,
           se,
           talpha0,
           tmu0,
           tSiRiSr0,
           tolerance,
           itermax,
           lnz_tol);
  }
  
  return(Rcpp::DataFrame::create(_["logodds"]=lovec,
                                 _["sigb"]=sigbvec,
                                 _["lnZ"]=nlzvec));
}

#endif






Rcpp::DataFrame grid_rss_varbvsr_serial(
    const c_Matrix_internal SiRiS,
    const c_arrayxd_internal sigma_beta,
    const c_arrayxd_internal logodds,
    const c_arrayxd_internal  betahat,
    const c_arrayxd_internal  se,
    const c_arrayxd_internal talpha0,
    const c_arrayxd_internal tmu0,
    const c_arrayxd_internal tSiRiSr0,
    double tolerance,
    int itermax,
    bool isVerbose,
    bool islnz_tol){
  
  using namespace Rcpp;
  size_t sigb_size= sigma_beta.size();
  size_t logodds_size=logodds.size();
  size_t tot_size=sigb_size*logodds_size;
  
  Rcpp::NumericVector nlzvec(tot_size);
  Rcpp::NumericVector sigbvec(tot_size);
  Rcpp::NumericVector lovec(tot_size);
  
  //  std::cout<<"Starting grid_rss_varbvsr (serial)"<<std::endl;
  
  Rcpp::LogicalVector verbose(1);
  verbose(0)=isVerbose;
  Rcpp::LogicalVector lnz_tol(1);
  lnz_tol(0)=islnz_tol;
  
  
  for(size_t t=0; t<tot_size; t++){
    size_t i=t%logodds_size;
    size_t j=t/logodds_size;
    sigbvec(t)=sigma_beta(j);
    lovec(t)=logodds(i);
    nlzvec(t)=rss_varbvsr_squarem_iter(SiRiS,
           sigma_beta(j),
           logodds(i),
           betahat,
           se,
           talpha0,
           tmu0,
           tSiRiSr0,
           tolerance,
           itermax,
           lnz_tol);
  }
  
  return(Rcpp::DataFrame::create(_["logodds"]=lovec,
                                 _["sigb"]=sigbvec,
                                 _["lnZ"]=nlzvec));
}










//[[Rcpp::export]]
double wrap_rss_varbvs_squarem_optim(Rcpp::NumericVector par, 
                                     const Matrix_external SiRiS,
                                     const arrayxd_external betahat,
                                     const arrayxd_external se,
                                     const arrayxd_external talpha0,
                                     const arrayxd_external tmu0,
                                     const arrayxd_external tSiRiSr0,
                                     double tolerance,
                                     int itermax,
                                     Rcpp::LogicalVector lnz_tol){
  return(-rss_varbvsr_squarem_iter(SiRiS,
                                   par(0),
                                   par(1),
                                   betahat,
                                   se,
                                   talpha0,
                                   tmu0,
                                   tSiRiSr0,
                                   tolerance,
                                   itermax,
                                   lnz_tol));
}




//[[Rcpp::export(name="rss_varbvs_squarem_trace")]]
Rcpp::NumericMatrix wrap_rss_varbvs_squarem_trace(
                                     const Matrix_external SiRiS,
                                     const double sigma_beta,
                                     const double logodds,
                                     const arrayxd_external betahat,
                                     const arrayxd_external se,
                                     const arrayxd_external talpha0,
                                     const arrayxd_external tmu0,
                                     const arrayxd_external tSiRiSr0,
                                     double tolerance,
                                     int itermax,
                                     Rcpp::LogicalVector lnz_tol){
  
  return(Rcpp::wrap(rss_varbvsr_squarem_iter_trace(SiRiS,sigma_beta,
                                                   logodds,betahat,se,talpha0,
                                                   tmu0,tSiRiSr0,tolerance,itermax,lnz_tol)));
}


//[[Rcpp::export]]
double wrap_rss_varbvs_naive_optim(Rcpp::NumericVector par, 
                                     const Matrix_external SiRiS,
                                     const arrayxd_external betahat,
                                     const arrayxd_external se,
                                     const arrayxd_external talpha0,
                                     const arrayxd_external tmu0,
                                     const arrayxd_external tSiRiSr0,
                                     double tolerance,
                                     int itermax,
                                     Rcpp::LogicalVector lnz_tol){
  

  return(-rss_varbvsr_naive_iter(SiRiS,
                                   par(0),
                                   par(1),
                                   betahat,
                                   se,
                                   talpha0,
                                   tmu0,
                                   tSiRiSr0,
                                   tolerance,
                                   itermax,
                                   Rcpp::LogicalVector::create(false),
                                   lnz_tol));
}







// 
// 
// #if RCPP_PARALLEL_USE_TBB
// 
// using namespace tbb;
// 
// 
// Rcpp::DataFrame fit_rss_logodds(
//     const c_sparseMatrix_internal SiRiS,
//     const c_arrayxd_internal sigma_beta,
//     const double logodds0,
//     const c_arrayxd_internal betahat,
//     const c_arrayxd_internal se,
//     const c_arrayxd_internal talpha0,
//     const c_arrayxd_internal tmu0,
//     const c_arrayxd_internal tSiRiSr0,
//     double tolerance,
//     int itermax,
//     Rcpp::LogicalVector verbose,
//     Rcpp::LogicalVector lnz_tol){
//   using namespace Rcpp;  
//   
//   size_t sigb_size= sigma_beta.size();
//   size_t tot_size=sigb_size;
//   Rcpp::NumericVector tlogodds(sigb_size);
//   std::fill(tlogodds.begin(),tlogodds.end(),logodds0);
//   Rcpp::NumericVector nlzvec(sigb_size);
//   //  Rcpp::NumericVector sigb_val(sigb_size);
//   Rcpp::NumericVector logodds_val(sigb_size);
//   
//   parallel_for(blocked_range<size_t>(0,tot_size),
//                [&](const blocked_range<size_t>& r){
//                  for(size_t t=r.begin(); t!=r.end(); t++){
//                    nlzvec[t]=rss_varbvsr_squarem_iter_fit_logodds(SiRiS,
//                                                       sigma_beta(t),
//                                                       tlogodds(t),
//                                                       betahat,
//                                                       se,
//                                                       talpha0,
//                                                       tmu0,
//                                                       tSiRiSr0,
//                                                       tolerance,
//                                                       itermax,
//                                                       lnz_tol);
//                    logodds_val(t)=tlogodds(t);
//                  }});
//   
//   
//   return DataFrame::create(_["logodds"]=logodds_val,
//                            _["sigb"]=sigma_beta,
//                            _["lnZ"]=nlzvec);
// }
// 
// 
// 
// 
// #else
// 
// Rcpp::DataFrame fit_rss_logodds(
//     const c_sparseMatrix_internal SiRiS,
//     const c_arrayxd_internal sigma_beta,
//     const double logodds0,
//     const c_arrayxd_internal betahat,
//     const c_arrayxd_internal se,
//     const c_arrayxd_internal talpha0,
//     const c_arrayxd_internal tmu0,
//     const c_arrayxd_internal tSiRiSr0,
//     double tolerance,
//     int itermax,
//     Rcpp::LogicalVector verbose,
//     Rcpp::LogicalVector lnz_tol){
//   using namespace Rcpp;  
//   
//   size_t sigb_size= sigma_beta.size();
//   size_t tot_size=sigb_size;
//   Rcpp::NumericVector tlogodds(sigb_size);
//   std::fill(tlogodds.begin(),tlogodds.end(),logodds0);
//   Rcpp::NumericVector nlzvec(sigb_size);
//   //  Rcpp::NumericVector sigb_val(sigb_size);
//   Rcpp::NumericVector logodds_val(sigb_size);
//   
//   for(size_t t=0; t<sigb_size; t++){
//     nlzvec[t]=rss_varbvsr_squarem_iter_fit_logodds(SiRiS,
//                                        sigma_beta(t),
//                                        tlogodds(t),
//                                        betahat,
//                                        se,
//                                        talpha0,
//                                        tmu0,
//                                        tSiRiSr0,
//                                        tolerance,
//                                        itermax,
//                                        lnz_tol);
//     logodds_val(t)=tlogodds(t);
//   }
//   
//   
//   return DataFrame::create(_["logodds"]=logodds_val,
//                            _["sigb"]=sigma_beta,
//                            _["lnZ"]=nlzvec);
// }
// 
// #endif
// 
// 
// 

// Rcpp::DataFrame rss_varbvsr_fit_hyperparameters(
//     const sparseMatrix_external SiRiS,
//     const arrayxd_external sigma_beta,
//     const arrayxd_external logodds0,
//     const arrayxd_external betahat,
//     const arrayxd_external se,
//     const arrayxd_external talpha0,
//     const arrayxd_external tmu0,
//     const arrayxd_external tSiRiSr0,
//     double tolerance,
//     int itermax,
//     Rcpp::LogicalVector verbose,
//     Rcpp::LogicalVector lnz_tol){
//   
//   return fit_rss_logodds(SiRiS,sigma_beta,logodds0(0),betahat,
//                           se,talpha0,tmu0,tSiRiSr0,tolerance,
//                           itermax,verbose,lnz_tol);
// }  

