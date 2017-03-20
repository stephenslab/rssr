#include <RcppEigen.h>
#include "rssr.h"
#include <cstdio>
#include <cmath>
//[[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>




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
  std::cout<<"Starting grid_rss_varbvsr (tbb)"<<std::endl;

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

#endif





#if RCPP_PARALLEL_USE_TBB

using namespace tbb;


Rcpp::DataFrame fit_rss_logodds(
    const c_sparseMatrix_internal SiRiS,
    const c_arrayxd_internal sigma_beta,
    const double logodds0,
    const c_arrayxd_internal betahat,
    const c_arrayxd_internal se,
    const c_arrayxd_internal talpha0,
    const c_arrayxd_internal tmu0,
    const c_arrayxd_internal tSiRiSr0,
    double tolerance,
    int itermax,
    Rcpp::LogicalVector verbose,
    Rcpp::LogicalVector lnz_tol){
  using namespace Rcpp;  
  
  size_t sigb_size= sigma_beta.size();
  size_t tot_size=sigb_size;
  Rcpp::NumericVector tlogodds(sigb_size);
  std::fill(tlogodds.begin(),tlogodds.end(),logodds0);
  Rcpp::NumericVector nlzvec(sigb_size);
  //  Rcpp::NumericVector sigb_val(sigb_size);
  Rcpp::NumericVector logodds_val(sigb_size);
  
  parallel_for(blocked_range<size_t>(0,tot_size),
               [&](const blocked_range<size_t>& r){
                 for(size_t t=r.begin(); t!=r.end(); t++){
                   nlzvec[t]=rss_varbvsr_squarem_iter_fit_logodds(SiRiS,
                                                      sigma_beta(t),
                                                      tlogodds(t),
                                                      betahat,
                                                      se,
                                                      talpha0,
                                                      tmu0,
                                                      tSiRiSr0,
                                                      tolerance,
                                                      itermax,
                                                      lnz_tol);
                   logodds_val(t)=tlogodds(t);
                 }});
  
  
  return DataFrame::create(_["logodds"]=logodds_val,
                           _["sigb"]=sigma_beta,
                           _["lnZ"]=nlzvec);
}




#else

Rcpp::DataFrame fit_rss_logodds(
    const c_sparseMatrix_internal SiRiS,
    const c_arrayxd_internal sigma_beta,
    const double logodds0,
    const c_arrayxd_internal betahat,
    const c_arrayxd_internal se,
    const c_arrayxd_internal talpha0,
    const c_arrayxd_internal tmu0,
    const c_arrayxd_internal tSiRiSr0,
    double tolerance,
    int itermax,
    Rcpp::LogicalVector verbose,
    Rcpp::LogicalVector lnz_tol){
  using namespace Rcpp;  
  
  size_t sigb_size= sigma_beta.size();
  size_t tot_size=sigb_size;
  Rcpp::NumericVector tlogodds(sigb_size);
  std::fill(tlogodds.begin(),tlogodds.end(),logodds0);
  Rcpp::NumericVector nlzvec(sigb_size);
  //  Rcpp::NumericVector sigb_val(sigb_size);
  Rcpp::NumericVector logodds_val(sigb_size);
  
  for(size_t t=0; t<sigb_size; t++){
    nlzvec[t]=rss_varbvsr_squarem_iter_fit_logodds(SiRiS,
                                       sigma_beta(t),
                                       tlogodds(t),
                                       betahat,
                                       se,
                                       talpha0,
                                       tmu0,
                                       tSiRiSr0,
                                       tolerance,
                                       itermax,
                                       lnz_tol);
    logodds_val(t)=tlogodds(t);
  }
  
  
  return DataFrame::create(_["logodds"]=logodds_val,
                           _["sigb"]=sigma_beta,
                           _["lnZ"]=nlzvec);
}

#endif



//[[Rcpp::export]]
Rcpp::DataFrame rss_varbvsr_fit_hyperparameters(
    const sparseMatrix_external SiRiS,
    const arrayxd_external sigma_beta,
    const arrayxd_external logodds0,
    const arrayxd_external betahat,
    const arrayxd_external se,
    const arrayxd_external talpha0,
    const arrayxd_external tmu0,
    const arrayxd_external tSiRiSr0,
    double tolerance,
    int itermax,
    Rcpp::LogicalVector verbose,
    Rcpp::LogicalVector lnz_tol){
  
  return fit_rss_logodds(SiRiS,sigma_beta,logodds0(0),betahat,
                          se,talpha0,tmu0,tSiRiSr0,tolerance,
                          itermax,verbose,lnz_tol);
}  

