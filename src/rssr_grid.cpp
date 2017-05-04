#include <RcppEigen.h>
#include "rssr.h"
#include <cstdio>
#include <cmath>
#include <RcppParallel.h>

#include <progress.hpp>
// [[Rcpp::depends(RcppProgress)]]
//[[Rcpp::depends(RcppParallel)]]


template<typename T> double rss_varbvsr_squarem_iter(const T SiRiS,
                                                     const double sigma_beta,
                                                     const double logodds,
                                                     const c_arrayxd_internal  betahat,
                                                     const c_arrayxd_internal se,
                                                     arrayxd_internal alpha,
                                                     arrayxd_internal mu,
                                                     arrayxd_internal SiRiSr,
                                                     const double tolerance,
                                                     const int itermax,
                                                     Rcpp::LogicalVector lnz_tol,int &iter, double &max_err){
  
  
  
  //This function implements RSS with variational bayes and the SQUAREM algorithm.
  iter=0;
  bool lnztol=lnz_tol[0];
  const size_t p = betahat.size();
  
  // Eigen::ArrayXd alpha=talpha0;
  // Eigen::ArrayXd mu=tmu0;
  // Eigen::ArrayXd SiRiSr=tSiRiSr0;
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
  
  max_err=1;
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












template<typename T> double rss_varbvsr_naive_iter(const T SiRiS,
                                                     const double sigma_beta,
                                                     const double logodds,
                                                     const c_arrayxd_internal  betahat,
                                                     const c_arrayxd_internal se,
                                                     arrayxd_internal alpha,
                                                     arrayxd_internal mu,
                                                     arrayxd_internal SiRiSr,
                                                     const double tolerance,
                                                     const int itermax,
                                                     Rcpp::LogicalVector lnz_tol,int &iter, double &max_err){
  
  
  
  
  //This function implements RSS with variational bayes and the SQUAREM algorithm.
  iter=0;
  bool lnztol=lnz_tol[0];
  const size_t p = betahat.size();
  
  // Eigen::ArrayXd alpha=talpha0;
  // Eigen::ArrayXd mu=tmu0;
  // Eigen::ArrayXd SiRiSr=tSiRiSr0;
  double lnZ=log(0);

  
  Eigen::ArrayXd talpha=alpha;
  Eigen::ArrayXd tmu=mu;
  Eigen::ArrayXd tSiRiSr=SiRiSr;
  
  double sigma_beta_square=sigma_beta*sigma_beta;
  Eigen::ArrayXd sesquare =se*se;
  Eigen::ArrayXd  q= betahat/sesquare;
  Eigen::ArrayXd  s= (sesquare*(sigma_beta*sigma_beta))/(sesquare+(sigma_beta*sigma_beta));
  Eigen::ArrayXd ssrat((s/sigma_beta_square).log());
  
  lnZ=calculate_lnZ(q,alpha*mu,SiRiSr,logodds,sesquare,alpha,mu,s,sigma_beta);
  
  
  max_err=1;
  double lnZ0=lnZ;
  double rel_l0=0;
  double rel_li=0;
  
  
  while(max_err>tolerance){
    talpha=alpha;
    tmu=mu;
    lnZ0=lnZ;
    bool reverse = iter%2!=0;
    rss_varbvsr_iter(SiRiS,sigma_beta_square,s,logodds,betahat,sesquare,ssrat,alpha,mu,SiRiSr,reverse);
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
    
    Eigen::ArrayXd alpha=talpha0;
    Eigen::ArrayXd mu=tmu0;
    Eigen::ArrayXd SiRiSr=tSiRiSr0;
    int iter=0;
    double max_err=1;
    
    double lnZ=rss_varbvsr_naive_iter<c_Matrix_internal>(SiRiS,sigma_beta,logodds,betahat,se,alpha,mu,SiRiSr,tolerance,itermax,lnz_tol,iter,max_err);
    
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
  
  Eigen::ArrayXd alpha=talpha0;
  Eigen::ArrayXd mu=tmu0;
  Eigen::ArrayXd SiRiSr=tSiRiSr0;
  int iter=0;
  double max_err=1;
  
  double lnZ=rss_varbvsr_naive_iter<c_sparseMatrix_internal>(SiRiS,sigma_beta,logodds,betahat,se,alpha,mu,SiRiSr,tolerance,itermax,lnz_tol,iter,max_err);
  
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
  
  
  Eigen::ArrayXd alpha=talpha0;
  Eigen::ArrayXd mu=tmu0;
  Eigen::ArrayXd SiRiSr=tSiRiSr0;
  
  int iter=0;
  double max_err=1;
  
  double lnZ=rss_varbvsr_squarem_iter<c_Matrix_internal>(SiRiS,sigma_beta,logodds,betahat,se,alpha,mu,SiRiSr,tolerance,itermax,lnz_tol,iter,max_err);
  
  
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
  
  
  Eigen::ArrayXd alpha=talpha0;
  Eigen::ArrayXd mu=tmu0;
  Eigen::ArrayXd SiRiSr=tSiRiSr0;
  int iter=0;
  double max_err=1;
  
  double lnZ=rss_varbvsr_squarem_iter<c_sparseMatrix_internal>(SiRiS,sigma_beta,logodds,betahat,se,alpha,mu,SiRiSr,tolerance,itermax,lnz_tol,iter,max_err);
  
  
  return Rcpp::List::create(Rcpp::Named("alpha")=alpha,
                            Rcpp::Named("mu")=mu,
                            Rcpp::Named("SiRiSr")=SiRiSr,
                            Rcpp::Named("max_err")=max_err,
                            Rcpp::Named("lnZ")=lnZ,
                            Rcpp::Named("iter")=iter);
}





#if RCPP_PARALLEL_USE_TBB

using namespace tbb;


template<typename T> Rcpp::DataFrame grid_rss_varbvsr(
    const T SiRiS,
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
  
  size_t tot_size=sigb_size;
  if(tot_size!=logodds_size){
    Rcpp::stop("Length of sigma_beta must equal length of logodds");
  }
  
  Eigen::ArrayXd npivec(tot_size);
  Eigen::ArrayXd nlzvec(tot_size);
  Eigen::ArrayXd sigbvec(tot_size);
  Eigen::ArrayXd lovec(tot_size);
  Eigen::ArrayXd pvevec(tot_size);
  
  
  Rcpp::LogicalVector verbose(1);
  verbose(0)=isVerbose;
  Rcpp::LogicalVector lnz_tol(1);
  lnz_tol(0)=islnz_tol;
  
  parallel_for(blocked_range<size_t>(0,tot_size),
               [&](const blocked_range<size_t>& r)  {
                 for(size_t t=r.begin(); t!=r.end(); t++){
                   size_t i=t;
                   size_t j=t;
                   sigbvec(t)=sigma_beta(j);
                   Eigen::ArrayXd copy_alpha(talpha0);
                   Eigen::ArrayXd copy_mu(tmu0);
                   Eigen::ArrayXd copy_SiRiSr(tSiRiSr0);
                   lovec(t)=logodds(i);
                   int iter=0;
                   double max_err=1;
                   double retvec=rss_varbvsr_squarem_iter(SiRiS,
                                                        sigma_beta(j),
                                                        logodds(i),
                                                        betahat,
                                                        se,
                                                        copy_alpha,
                                                        copy_mu,
                                                        copy_SiRiSr,
                                                        tolerance,
                                                        itermax,
                                                        lnz_tol,iter,max_err);
                   nlzvec[t]=retvec;
                   npivec[t]=copy_alpha.mean();
                   pvevec[t]=copy_SiRiSr.matrix().transpose()*(copy_alpha*copy_mu).matrix();
                 }});
  
  //  Rcpp::Rcout<<"mean lnZ is: "<<mean(nlzvec)<<std::endl;
  return(Rcpp::DataFrame::create(_["logodds"]=lovec,
                                 _["sigb"]=sigbvec,
                                 _["alpha_mean"]=npivec,
                                 _["pve"]=pvevec,
                                 _["lnZ"]=nlzvec));
}




#else


templtate<typename T> Rcpp::DataFrame grid_rss_varbvsr(
    const T SiRiS,
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
  size_t tot_size=sigb_size;
  if(tot_size!=logodds_size){
    Rcpp::stop("Length of sigma_beta must equal length of logodds");
  }
  
  Eigen::ArrayXd npivec(tot_size);
  Eigen::ArrayXd pvevec(tot_size);
  Rcpp::NumericVector nlzvec(tot_size);
  Rcpp::NumericVector sigbvec(tot_size);
  Rcpp::NumericVector lovec(tot_size);
  
  //   std::cout<<"Starting grid_rss_varbvsr (serial)"<<std::endl;
  
  Rcpp::LogicalVector verbose(1);
  verbose(0)=isVerbose;
  Rcpp::LogicalVector lnz_tol(1);
  lnz_tol(0)=islnz_tol;
  
  
  for(size_t t=0; t<tot_size; t++){
    Eigen::ArrayXd copy_alpha(talpha0);
    Eigen::ArrayXd copy_mu(tmu0);
    Eigen::ArrayXd copy_SiRiSr(tSiRiSr0);
    size_t i=t;
    size_t j=t;
    sigbvec(t)=sigma_beta(j);
    lovec(t)=logodds(i);
    int iter=0;
    double max_err=1;
    
    
    double retvec =rss_varbvsr_squarem_iter(SiRiS,
                                          sigma_beta(j),
                                          logodds(i),
                                          betahat,
                                          se,
                                          copy_alpha,
                                          copy_mu,
                                          copy_SiRiSr,
                                          tolerance,
                                          itermax,
                                          lnz_tol,iter,max_err);
    nlzvec[t]=retvec[0];
    npivec[t]=copy_alpha.mean();
    pvevec[t]=copy_SiRiSr.matrix().transpose()*(copy_alpha*copy_mu).matrix();
  }
  return(Rcpp::DataFrame::create(_["logodds"]=lovec,
                                 _["sigb"]=sigbvec,
                                 _["alpha_mean"]=npivec,
                                 _["pve"]=pvevec,
                                 _["lnZ"]=nlzvec));
}

#endif



//[[Rcpp::export]]
Rcpp::DataFrame grid_search_rss_varbvsr_sp(
    const  sparseMatrix_external SiRiS,
    const arrayxd_external sigma_beta,
    const arrayxd_external logodds,
    const arrayxd_external  betahat,
    const arrayxd_external  se,
    const arrayxd_external talpha0,
    const arrayxd_external tmu0,
    const arrayxd_external tSiRiSr0,
    double tolerance,
    int itermax,
    Rcpp::LogicalVector verbose,
    Rcpp::LogicalVector lnz_tol){
  //  std::cout<<"Starting grid_search_rss_varbvsr"<<std::endl;
  bool isVerbose = verbose(0);
  bool islnz_tol = lnz_tol(0);
  
  return grid_rss_varbvsr<c_sparseMatrix_internal>(SiRiS,sigma_beta,logodds,betahat,
                          se,talpha0,tmu0,tSiRiSr0,tolerance,
                          itermax,isVerbose,islnz_tol);
} 


//[[Rcpp::export]]
Rcpp::DataFrame grid_search_rss_varbvsr(
    const  Matrix_external SiRiS,
    const arrayxd_external sigma_beta,
    const arrayxd_external logodds,
    const arrayxd_external  betahat,
    const arrayxd_external  se,
    const arrayxd_external talpha0,
    const arrayxd_external tmu0,
    const arrayxd_external tSiRiSr0,
    double tolerance,
    int itermax,
    Rcpp::LogicalVector verbose,
    Rcpp::LogicalVector lnz_tol){
  //  std::cout<<"Starting grid_search_rss_varbvsr"<<std::endl;
  bool isVerbose = verbose(0);
  bool islnz_tol = lnz_tol(0);
  
  if(sigma_beta.minCoeff()<=0){
    Rcpp::stop("sigma_beta must be strictly positive");
  }
  
  return grid_rss_varbvsr<c_Matrix_internal>(SiRiS,sigma_beta,logodds,betahat,
                          se,talpha0,tmu0,tSiRiSr0,tolerance,
                          itermax,isVerbose,islnz_tol);
} 






