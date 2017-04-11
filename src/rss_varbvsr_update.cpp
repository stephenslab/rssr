// For a description of this C code, see rss_varbvsr_update.m.
#include <RcppEigen.h>
#include "rssr.h"

#include <cstdio>





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
  
  
  rss_varbvsr_iter(SiRiS,
                   sigma_beta,
                   logodds,
                   betahat,
                   se,
                   talpha,
                   tmu,
                   tSiRiSr,
                   reverse);
  return Rcpp::List::create(Rcpp::Named("alpha1")=talpha,
                            Rcpp::Named("mu1")=tmu,
                            Rcpp::Named("SiRiSr")=tSiRiSr);
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
  
  
  rss_varbvsr_iter(SiRiS,
                   sigma_beta,
                   logodds,
                   betahat,
                   se,
                   talpha,
                   tmu,
                   tSiRiSr,
                   reverse);
  return Rcpp::List::create(Rcpp::Named("alpha1")=talpha,
                            Rcpp::Named("mu1")=tmu,
                            Rcpp::Named("SiRiSr")=tSiRiSr);
}











//' Run RSS with the variational bayes algorithm accelerated with SQUAREM
//' @template rssr
//' @param talpha0 a length p vector specifying the initial value of alpha
//' @param tmu0 a length p vector specifying the initial value of mu
//' @param SiRiSr0 a length p vector specifying the initial value of SiRiSr
//[[Rcpp::export]]
Rcpp::List rss_varbvsr_squarem_sp(const sparseMatrix_external SiRiS,
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
  
  
  //This function implements RSS with variational bayes and the SQUAREM algorithm.
  
  bool verbosev=verbose[0];
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
      // Rcpp::Rcout<<"mtp <1: "<<mtp<<std::endl;
      // Rcpp::Rcout<<"mean(alpha.col(c)): "<<alpha.mean()<<std::endl;
      // Rcpp::Rcout<<"mean(mu.col(c)): "<<mu.mean()<<std::endl;
      // Rcpp::Rcout<<"mean(SiRiSr.col(c)): "<<SiRiSr.mean()<<std::endl;
      alpha=alpha0-2*mtp*alpha_r+(mtp*mtp)*alpha_v;
      mu=mu0-2*mtp*mu_r+(mtp*mtp)*mu_v;
      SiRiSr=SiRiS*(alpha*mu).matrix();
      // Rcpp::Rcout<<"mean(alpha.col(c)): "<<alpha.mean()<<std::endl;
      // Rcpp::Rcout<<"mean(mu.col(c)): "<<mu.mean()<<std::endl;
      // Rcpp::Rcout<<"mean(SiRiSr.col(c)): "<<SiRiSr.mean()<<std::endl;
    }
    
    rss_varbvsr_iter(SiRiS,sigma_beta,logodds,betahat,se,alpha,mu,SiRiSr,reverse);
    lnZ=  calculate_lnZ(q,alpha*mu,SiRiSr,logodds,sesquare,alpha,mu,s,sigma_beta);
    if((mtp<(-1)) && (lnZ < lnZ0)){
      Rcpp::Rcout<<"begin bt "<<std::endl;
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
    if(verbosev){
      double absr=(alpha*mu).abs().maxCoeff();
      int asum=round(alpha.sum());
      printf("%4d %+13.6e %1.9e %4d %1.9e %1.9e\n",(int)iter,lnZ,max_err,(int) asum,rel_l0,rel_li);
    }
    if(iter>itermax){
      printf("Maximum iteration number reached: %+0.2d \n",(int)iter);
      printf("The log variational lower bound of the last step increased by %+0.2e\n",lnZ-lnZ0);
      break;
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


//' Run RSS with the variational bayes algorithm accelerated with SQUAREM
//' @template rssr
//' @param talpha0 a length p vector specifying the initial value of alpha
//' @param tmu0 a length p vector specifying the initial value of mu
//' @param SiRiSr0 a length p vector specifying the initial value of SiRiSr
//[[Rcpp::export]]
Rcpp::List rss_varbvsr_squarem(const Matrix_external SiRiS,
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
  
  
  //This function implements RSS with variational bayes and the SQUAREM algorithm.
  
  bool verbosev=verbose[0];
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
      // Rcpp::Rcout<<"mtp <1: "<<mtp<<std::endl;
      // Rcpp::Rcout<<"mean(alpha.col(c)): "<<alpha.mean()<<std::endl;
      // Rcpp::Rcout<<"mean(mu.col(c)): "<<mu.mean()<<std::endl;
      // Rcpp::Rcout<<"mean(SiRiSr.col(c)): "<<SiRiSr.mean()<<std::endl;
      alpha=alpha0-2*mtp*alpha_r+(mtp*mtp)*alpha_v;
      mu=mu0-2*mtp*mu_r+(mtp*mtp)*mu_v;
      SiRiSr=SiRiS*(alpha*mu).matrix();
      // Rcpp::Rcout<<"mean(alpha.col(c)): "<<alpha.mean()<<std::endl;
      // Rcpp::Rcout<<"mean(mu.col(c)): "<<mu.mean()<<std::endl;
      // Rcpp::Rcout<<"mean(SiRiSr.col(c)): "<<SiRiSr.mean()<<std::endl;
    }
    
    rss_varbvsr_iter(SiRiS,sigma_beta,logodds,betahat,se,alpha,mu,SiRiSr,reverse);
    lnZ=  calculate_lnZ(q,alpha*mu,SiRiSr,logodds,sesquare,alpha,mu,s,sigma_beta);
    if((mtp<(-1)) && (lnZ < lnZ0)){
      Rcpp::Rcout<<"begin bt "<<std::endl;
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
    if(verbosev){
      double absr=(alpha*mu).abs().maxCoeff();
      int asum=round(alpha.sum());
      printf("%4d %+13.6e %1.9e %4d %1.9e %1.9e\n",(int)iter,lnZ,max_err,(int) asum,rel_l0,rel_li);
    }
    if(iter>itermax){
      printf("Maximum iteration number reached: %+0.2d \n",(int)iter);
      printf("The log variational lower bound of the last step increased by %+0.2e\n",lnZ-lnZ0);
      break;
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




//' Estimate hyperparameters of RSS using the EM algorithm
//' @template rssr
//' @param talpha0 a length p vector specifying the initial value of alpha
//' @param tmu0 a length p vector specifying the initial value of mu
//' @param SiRiSr0 a length p vector specifying the initial value of SiRiSr
//[[Rcpp::export]]
Rcpp::List rss_varbvsr_squarem_fit_logodds(const sparseMatrix_external SiRiS,
                                           const double sigma_beta,
                                           const double logodds0,
                                           const arrayxd_external betahat,
                                           const arrayxd_external se,
                                           const arrayxd_external talpha0,
                                           const arrayxd_external tmu0,
                                           const arrayxd_external tSiRiSr0,
                                           double tolerance,
                                           int itermax,
                                           Rcpp::LogicalVector verbose,
                                           Rcpp::LogicalVector lnz_tol){
  
  
  //This function implements RSS with variational bayes and the SQUAREM algorithm.
  
  bool verbosev=verbose[0];
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
  
  double logodds=logodds0;
  double pi=0;
  
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
    if(verbosev){
      double absr=(alpha*mu).abs().maxCoeff();
      int asum=round(alpha.sum());
      printf("%4d %+13.6e %1.9e %4d %1.9e %1.9e\n",(int)iter,lnZ,max_err,(int) asum,rel_l0,rel_li);
    }
    
    logodds=update_logodds(alpha);

    if(iter>itermax){
      printf("Maximum iteration number reached: %+0.2d \n",(int)iter);
      printf("The log variational lower bound of the last step increased by %+0.2e\n",lnZ-lnZ0);
      break;
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
Rcpp::List rss_varbvsr_naive_sp (const sparseMatrix_external SiRiS,
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
  
  return Rcpp::List::create(Rcpp::Named("alpha")=alpha,
                            Rcpp::Named("mu")=mu,
                            Rcpp::Named("SiRiSr")=SiRiSr,
                            Rcpp::Named("max_err")=max_err,
                            Rcpp::Named("lnZ")=lnZ,
                            Rcpp::Named("iter")=iter);
  
}




//[[Rcpp::export]]
Rcpp::List rss_varbvsr_naive (const Matrix_external SiRiS,
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
  
  return Rcpp::List::create(Rcpp::Named("alpha")=alpha,
                            Rcpp::Named("mu")=mu,
                            Rcpp::Named("SiRiSr")=SiRiSr,
                            Rcpp::Named("max_err")=max_err,
                            Rcpp::Named("lnZ")=lnZ,
                            Rcpp::Named("iter")=iter);
  
}
