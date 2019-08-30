#include <RcppEigen.h>
#include "rssr.h"
#include <cstdio>
#include <cmath>
#include <RcppParallel.h>
//[[Rcpp::depends(RcppParallel)]]






template <typename T> void squarem_backtrack(arrayxd_internal &alpha,
                                             const c_arrayxd_internal alpha0,
                                             const c_arrayxd_internal alpha1,
                                             arrayxd_internal &mu,
                                             const c_arrayxd_internal mu0,
                                             const c_arrayxd_internal mu1,
                                             arrayxd_internal &SiRiSr,
                                             T SiRiS,
                                             const double sigma_beta,
                                             const double logodds,
                                             const c_arrayxd_internal s,
                                             const c_arrayxd_internal betahat,
                                             const c_arrayxd_internal sesquare,
                                             const c_arrayxd_internal ssrat,
                                             double &mtp, double &lnZ, double &lnZ0,double &lnZ00, bool &doDie)
{
  
  
  // alpha_r=alpha1-alpha0;
  // alpha_v=(alpha-alpha1)-alpha_r;
  // 
  // mu_r   = mu1-mu0;
  // mu_v   = mu-mu1-mu_r;
//  Rcpp::Rcerr<<"Do we need to backtrack?"<<std::endl;
  double sigma_beta_square=sigma_beta*sigma_beta;
  
  size_t max_backtrack=12;
  if((mtp<(-1)) && (lnZ < lnZ0)){
    size_t num_bt=0;
    // Rcpp::Rcerr<<"Starting backtracking"<<std::endl;
    while((lnZ<lnZ0) &&(num_bt < max_backtrack)){
      mtp = 0.5*(mtp-1);
      alpha = alpha0-2*mtp*(alpha1-alpha0)+(mtp*mtp)*(alpha-2*alpha1+alpha0);
      mu = mu0-2*mtp*(mu1-mu0)+(mtp*mtp)*(mu-2*mu1+mu0);
      SiRiSr = (SiRiS*(alpha*mu).matrix()).array();
      rss_varbvsr_iter(SiRiS,sigma_beta_square,s,logodds,betahat,sesquare,ssrat,alpha,mu,SiRiSr,false);
      rss_varbvsr_iter(SiRiS,sigma_beta_square,s,logodds,betahat,sesquare,ssrat,alpha,mu,SiRiSr,true);
      // r=alpha*mu;
      lnZ=calculate_lnZ(betahat/sesquare,alpha*mu,SiRiSr,logodds,sesquare,alpha,mu,s,sigma_beta);          
      num_bt=num_bt+1;
      if(num_bt==max_backtrack){
        alpha=alpha0;
        mu=mu0;
        SiRiSr = (SiRiS*(alpha*mu).matrix()).array();
        lnZ0=calculate_lnZ(betahat/sesquare,alpha*mu,SiRiSr,logodds,sesquare,alpha,mu,s,sigma_beta);
        rss_varbvsr_iter(SiRiS,sigma_beta_square,s,logodds,betahat,sesquare,ssrat,alpha,mu,SiRiSr,false);
        rss_varbvsr_iter(SiRiS,sigma_beta_square,s,logodds,betahat,sesquare,ssrat,alpha,mu,SiRiSr,true);
        lnZ=calculate_lnZ(betahat/sesquare,alpha*mu,SiRiSr,logodds,sesquare,alpha,mu,s,sigma_beta);
        if(lnZ<lnZ0){
          doDie=true;
        }
        mtp=-1;
      }
      // Rcpp::Rcerr<<"Final!: num_bt:"<<num_bt<<" step_size:"<<mtp<<"lnZ:"<<lnZ<<" lnZ0:"<<lnZ0<<std::endl;
    }
  }

}


double calculate_mtp(const c_arrayxd_internal alpha,const c_arrayxd_internal alpha1,const c_arrayxd_internal alpha0, const c_arrayxd_internal mu,const c_arrayxd_internal mu1, const c_arrayxd_internal mu0){
  return(  -std::sqrt(((alpha1-alpha0).square()).sum()+((mu1-mu0).square()).sum())/std::sqrt(((alpha-2*alpha1+alpha0).square()).sum()+((mu-2*mu1+mu0).square()).sum())+double_lim::epsilon());
}

double calculate_mtp(const c_arrayxd_internal alpha_r,const c_arrayxd_internal alpha_v,const c_arrayxd_internal mu_r, const c_arrayxd_internal mu_v){
  
  return(  -std::sqrt(((alpha_r).square()).sum()+((mu_r).square()).sum())/std::sqrt(((alpha_v).square()).sum()+((mu_v).square()).sum())+double_lim::epsilon());
}

template <typename T> void squarem_adjust(arrayxd_internal alpha,
                                          const c_arrayxd_internal alpha0,
                                          const c_arrayxd_internal alpha1,
                                          const c_arrayxd_internal alpha_r,
                                          const c_arrayxd_internal alpha_v,
                                          arrayxd_internal mu,
                                          const c_arrayxd_internal mu0,
                                          const c_arrayxd_internal mu1,
                                          const c_arrayxd_internal mu_r,
                                          const c_arrayxd_internal mu_v,
                                          arrayxd_internal SiRiSr,
                                          T SiRiS,
                                          double &mtp){
  

  
  mtp=calculate_mtp(alpha_r,alpha_v,mu_r,mu_v);
 
  if(mtp > (-1)){
    mtp=-1;
  }else{
   if(!std::isnan(mtp)){
    alpha=alpha0-2*mtp*(alpha_r)+(mtp*mtp)*(alpha_v);
    mu=mu0-2*mtp*(mu_r)+(mtp*mtp)*(mu_v);
    SiRiSr=(SiRiS*(alpha*mu).matrix()).array();
  }
  }
  
}







void calc_max_err(const double lnZ,const double lnZ0,const c_arrayxd_internal alpha, const c_arrayxd_internal alpha0,const c_arrayxd_internal mu,const c_arrayxd_internal mu0,  double &max_err,const bool lnztol){
  
  if(lnztol){
    max_err=rel_err(lnZ,lnZ0);
  }else{
    max_err=find_maxerr(alpha,alpha0,alpha*mu,alpha0*mu0);
  }
}

double initialize_value(const double val, const double param){
  return(val);
}



Eigen::ArrayXd initialize_array(const c_arrayxd_internal init, const double param){
  Eigen::ArrayXd ret=init;
  return(ret);
}



Eigen::ArrayXd initialize_s(const c_arrayxd_internal sesquare, const double sigma_beta_square){
  return((sesquare*(sigma_beta_square))/(sesquare+(sigma_beta_square)));
}



Eigen::ArrayXd initialize_ssrat(const c_arrayxd_internal s, const double sigma_beta_square){
  return((s/sigma_beta_square).log());
}



//[[Rcpp::export]]
Rcpp::List wrap_squarem_adjust_prep(const Matrix_external SiRiS,
                                    const double sigma_beta,
                                    const double logodds,
                                    const arrayxd_external  betahat,
                                    const arrayxd_external se,
                                    arrayxd_external talpha,
                                    arrayxd_external tmu,
                                    arrayxd_external tSiRiSr,
                                    const double tolerance,
                                    const int itermax,
                                    Rcpp::LogicalVector lnz_tol){
  

  
  Eigen::ArrayXd alpha=talpha;
  Eigen::ArrayXd mu=tmu;
  Eigen::ArrayXd SiRiSr=tSiRiSr;
  
  Eigen::ArrayXd alpha2=talpha;
  Eigen::ArrayXd mu2=tmu;
  Eigen::ArrayXd SiRiSr2=tSiRiSr;
  
  double max_err=1;
  int iter=0;
  bool lnztol=lnz_tol[0];
  const size_t p = betahat.size();

  double lnZ= 0;
  
  
  Eigen::ArrayXd alpha0 = alpha;
  Eigen::ArrayXd mu0 =mu;
  Eigen::ArrayXd  SiRiSr0 = SiRiSr;

  ArrayXd alpha1=alpha;
  ArrayXd mu1=mu;
  ArrayXd SiRiSr1=SiRiSr;

  ArrayXd r=alpha*mu;
  
  
  double sigma_beta_square=sigma_beta*sigma_beta;
  Eigen::ArrayXd sesquare =se*se;
  Eigen::ArrayXd  q= betahat/sesquare;
  ArrayXd s=initialize_s(sesquare,sigma_beta_square);
  ArrayXd ssrat=initialize_ssrat(s,sigma_beta_square);

  if(r.hasNaN()){
    Rcpp::Rcerr<<"In iteration iter 0(0)"<<std::endl;
    Rcpp::stop("alpha*mu is not finite!");
  }
  
  lnZ=calculate_lnZ(q,r,SiRiSr,logodds,sesquare,alpha,mu,s,sigma_beta);
  
  double mtp=0;
  
  max_err=1;
  double lnZ0=lnZ;

  lnZ0=lnZ;
  alpha0=alpha;
  mu0=mu; 
  
  bool reverse = iter%2!=0;
  rss_varbvsr_iter(SiRiS,sigma_beta_square,s,logodds,betahat,sesquare,ssrat,alpha,mu,SiRiSr,reverse);
  
  alpha1=alpha;
  mu1=mu;
  SiRiSr1=SiRiSr;
  
  rss_varbvsr_iter(SiRiS,sigma_beta_square,s,logodds,betahat,sesquare,ssrat,alpha,mu,SiRiSr,reverse);

  alpha2=alpha;
  mu2=mu;
  SiRiSr2=SiRiSr;
  ArrayXd alpha_r=(alpha1-alpha0);
  ArrayXd alpha_v=(alpha-alpha1)-alpha_r;
  
  ArrayXd mu_r=(mu1-mu0);
  ArrayXd mu_v=(mu-mu1)-mu_r;
  squarem_adjust(alpha,alpha0,alpha1,alpha_r,alpha_v,mu,mu0,mu1,mu_r,mu_v,tSiRiSr,SiRiS,mtp);
  
  mtp= -std::sqrt(((alpha1-alpha0).square()).sum()+((mu1-mu0).square()).sum())/std::sqrt(((alpha-2*alpha1+alpha0).square()).sum()+((mu-2*mu1+mu0).square()).sum());
  
  // Rcpp::Rcout<<"squarem adjust"<<std::endl;
  alpha=alpha0-2*mtp*(alpha1-alpha0)+(mtp*mtp)*(alpha-2*alpha1+alpha0);
  mu=mu0-2*mtp*(mu1-mu0)+(mtp*mtp)*(mu-2*mu1+mu0);
  SiRiSr=SiRiS*(alpha*mu).matrix();
  
  
  
  return(Rcpp::List::create(Rcpp::Named("alpha0")=alpha0,
                            Rcpp::Named("alpha1")=alpha1,
                            Rcpp::Named("alpha2")=alpha2,
                            Rcpp::Named("alpha")=alpha,
                            Rcpp::Named("mu0")=mu0,
                            Rcpp::Named("mu1")=mu1,
                            Rcpp::Named("mu2")=mu2,
                            Rcpp::Named("mu")=mu,
                            Rcpp::Named("SiRiSr")=tSiRiSr,
                            Rcpp::Named("mtp")=mtp));
}


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
  
  using namespace Eigen;

  //This function implements RSS with variational bayes and the SQUAREM algorithm.
  iter=0;
  bool lnztol=lnz_tol[0];
  const size_t p = betahat.size();
  

  double lnZ= initialize_value(0,logodds);

  
  
  Eigen::ArrayXd alpha0 = alpha;
  Eigen::ArrayXd mu0 =mu;
  Eigen::ArrayXd  SiRiSr0 = SiRiSr;
  
  
  ArrayXd alpha1=alpha;
  ArrayXd mu1=mu;
  ArrayXd SiRiSr1=SiRiSr;
  
  
  ArrayXd r=alpha*mu;
  
  
  double sigma_beta_square=sigma_beta*sigma_beta;
  Eigen::ArrayXd sesquare =se*se;
  Eigen::ArrayXd  q= betahat/sesquare;
  ArrayXd s=initialize_s(sesquare,sigma_beta_square);
  ArrayXd ssrat=initialize_ssrat(s,sigma_beta_square);

  
  ArrayXd alpha_r=(alpha1-alpha0);
  ArrayXd alpha_v=(alpha-alpha1)-alpha_r;
  
  ArrayXd mu_r=(mu1-mu0);
  ArrayXd mu_v=(mu-mu1)-mu_r;
  if(r.hasNaN()){
    Rcpp::Rcerr<<"In iteration iter 0(0)"<<std::endl;
    Rcpp::stop("alpha*mu is not finite!");
  }
  
  lnZ=calculate_lnZ(q,r,SiRiSr,logodds,sesquare,alpha,mu,s,sigma_beta);
  
  double mtp=-1;
  
  max_err=1;
  double lnZ0=lnZ;
  double lnZ00=lnZ0;
  
  while(max_err>tolerance){
    lnZ00=lnZ0;
    lnZ0=lnZ;
    alpha0=alpha;
    mu0=mu; 
    
    bool reverse = iter%2!=0;
    
    rss_varbvsr_iter(SiRiS,sigma_beta_square,s,logodds,betahat,sesquare,ssrat,alpha,mu,SiRiSr,reverse);
    
    
    alpha1=alpha;
    mu1=mu;
    SiRiSr1=SiRiSr;

    rss_varbvsr_iter(SiRiS,sigma_beta_square,s,logodds,betahat,sesquare,ssrat,alpha,mu,SiRiSr,reverse);

    alpha_r=(alpha1-alpha0);
    alpha_v=(alpha-alpha1)-alpha_r;
    
    mu_r=(mu1-mu0);
    mu_v=(mu-mu1)-mu_r;
    squarem_adjust<T>(alpha,alpha0,alpha1,alpha_r,alpha_v,mu,mu0,mu1,mu_r,mu_v,SiRiSr,SiRiS,mtp);

    rss_varbvsr_iter(SiRiS,sigma_beta_square,s,logodds,betahat,sesquare,ssrat,alpha,mu,SiRiSr,reverse);

    lnZ=  calculate_lnZ(q,alpha*mu,SiRiSr,logodds,sesquare,alpha,mu,s,sigma_beta);
    bool doDie=false;
    squarem_backtrack<T>(alpha,alpha0,alpha1,
                         mu,mu0,mu1,
                         SiRiSr,SiRiS,
                         sigma_beta,logodds,
                         s,betahat,
                         sesquare,ssrat,
                         mtp,lnZ, lnZ0,lnZ00, doDie);
    if(doDie){
      iter=itermax+1;
    }
    
    calc_max_err(lnZ,lnZ0, alpha, alpha0, mu, mu0,max_err,lnztol);
    
    if(iter>itermax){
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
  double lnZ=1;
  
  
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
  return(lnZ);
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
    const T tSiRiS,
    const c_arrayxd_internal tsigma_beta,
    const c_arrayxd_internal tlogodds,
    const c_arrayxd_internal  tbetahat,
    const c_arrayxd_internal  tse,
    const c_arrayxd_internal talpha0,
    const c_arrayxd_internal tmu0,
    const c_arrayxd_internal tSiRiSr0,
    double tolerance,
    int itermax,
    bool isVerbose,
    bool islnz_tol){
  
  using namespace Rcpp;
  using namespace Eigen;
  
  ArrayXd sigma_beta=tsigma_beta;
  ArrayXd logodds=tlogodds;
  ArrayXd betahat=tbetahat;
  ArrayXd se=tse;
  MatrixXd SiRiS=tSiRiS;

  
  
  
  size_t sigb_size= sigma_beta.size();
  size_t logodds_size=logodds.size();
  
  size_t tot_size=sigb_size;
  if(tot_size!=logodds_size){
    Rcpp::stop("Length of sigma_beta must equal length of logodds");
  }

  Eigen::ArrayXd npivec(tot_size);
  Eigen::ArrayXd nlzvec(tot_size);
  Eigen::ArrayXd errvec(tot_size);
  Eigen::ArrayXi itervec(tot_size);
  Eigen::ArrayXd pvevec(tot_size);
  size_t p=betahat.size();
  
  Rcpp::LogicalVector verbose(1);
  verbose(0)=isVerbose;
  Rcpp::LogicalVector lnz_tol(1);
  lnz_tol(0)=islnz_tol;
  // static affinity_partitioner ap;
  Eigen::ArrayXXd  all_alphas(p,tot_size);
  Eigen::ArrayXXd  all_mus(p,tot_size);
  Eigen::ArrayXXd  all_SiRiSr(p,tot_size);
  
  for(size_t i=0; i<tot_size; i++){
    all_alphas.col(i)=talpha0;
    all_mus.col(i)=tmu0;
    all_SiRiSr.col(i)=tSiRiSr0;
  }
  
  
  parallel_for(blocked_range<size_t>(0,tot_size),
               [&](const blocked_range<size_t>& r)  {
                 for(size_t t=r.begin(); t!=r.end(); t++){
                   size_t i=t;
                   size_t j=t;
                   int iter=0;
                   double max_err=1;
                   double retvec=rss_varbvsr_squarem_iter(SiRiS,
                                                          sigma_beta(j),
                                                          logodds(i),
                                                          betahat,
                                                          se,
                                                          all_alphas.col(t),
                                                          all_mus.col(t),
                                                          all_SiRiSr.col(t),
                                                          tolerance,
                                                          itermax,
                                                          lnz_tol,iter,max_err);
                   nlzvec[t]=retvec;
                   itervec[t]=iter;
                   errvec[t]=max_err;
                   npivec[t]=all_alphas.col(t).mean();
                   pvevec[t]=all_SiRiSr.col(t).matrix().transpose()*(all_alphas.col(t)*all_alphas.col(t)).matrix();
                 }});
  
  //  Rcpp::Rcout<<"mean lnZ is: "<<mean(nlzvec)<<std::endl;
  return(Rcpp::DataFrame::create(_["logodds"]=logodds,
                                 _["sigb"]=sigma_beta,
                                 _["rel_err"]=errvec,
                                 _["iterations"]=itervec,
                                 _["alpha_mean"]=npivec,
                                 _["pve"]=pvevec,
                                 _["lnZ"]=nlzvec));
}


// 
// Rcpp::DataFrame grid_rss_varbvsr_multitrait(
//     const Rcpp::NumericMatrix &R,
//     const c_arrayxd_internal sigma_beta,
//     const c_arrayxd_internal logodds,
//     const Rcpp::NumericMatrix  &betahat,
//     const Rcpp::NumericMatrix  &se,
//     const c_arrayxd_internal talpha0,
//     const c_arrayxd_internal tmu0,
//     double tolerance,
//     int itermax,
//     bool isVerbose,
//     const c_arrayxi_internal fgeneid,
//     bool islnz_tol){
//   
//   using namespace Rcpp;
//   size_t ngenes=betahat.cols();
//   size_t sigb_size= sigma_beta.size();
//   size_t logodds_size=logodds.size();
//   
//   if(ngenes!=se.cols()){
//     Rcpp::stop("se and betahat must have same dimensions!");
//   }
//   
//   size_t tot_size=sigb_size;
//   if(tot_size!=logodds_size){
//     Rcpp::stop("Length of sigma_beta must equal length of logodds");
//   }
//   
//   size_t ntot_size=tot_size*ngenes;
//   Eigen::ArrayXd npivec(ngenes*tot_size);
//   Eigen::ArrayXd nlzvec(ngenes*tot_size);
//   Eigen::ArrayXd errvec(ngenes*tot_size);
//   Eigen::ArrayXi itervec(ngenes*tot_size);
//   Eigen::ArrayXd pvevec(ngenes*tot_size);
//   Eigen::ArrayXi fgeneidmat(ngenes*tot_size);
//   
//   Eigen::ArrayXd sigbvec(ngenes*tot_size);
//   Eigen::ArrayXd logoddsvec(ngenes*tot_size);
// 
//   
//   
//   Rcpp::LogicalVector verbose(1);
//   verbose(0)=isVerbose;
//   Rcpp::LogicalVector lnz_tol(1);
//   lnz_tol(0)=islnz_tol;
//   // static affinity_partitioner ap;
//   parallel_for(blocked_range<size_t>(0,ngenes),
//                [&](const blocked_range<size_t>& r)  {
//                  for( size_t g=r.begin(); g!=r.end(); g++){
//                    // parallel_for(blocked_range<size_t>(0,tot_size),
//                    // for(size_t g=0;g<ngenes;g++){
//                    
//                    Eigen::VectorXd si = 1/se.col(g).array();
//                    Eigen::MatrixXd siris = si.asDiagonal()*R*si.asDiagonal();
//                    Eigen::ArrayXd sirisr = siris*(talpha0*tmu0).matrix();                    
//                    // for(size_t t=0;t<tot_size;t++){
//                    for(size_t t=0; t<tot_size; t++){
//                      // std::cout<<"t:"<<t<<std::endl;
//                      //   std::cout<<"g:"<<g<<std::endl;
//                      size_t i=t;
//                      size_t j=t;
//                      size_t mat_ind=(ngenes)*(g)+t;
//                      
//                      // Rcpp::Rcout<<"mat_ind:"<<mat_ind<<std::endl;
//                      
//                      Eigen::ArrayXd copy_alpha(talpha0);
//                      Eigen::ArrayXd copy_mu(tmu0);
//                      Eigen::ArrayXd copy_SiRiSr(sirisr);
//                      //                   lovec(t)=logodds(i);
//                      int iter=0;
//                      double max_err=1;
//                      double retvec=rss_varbvsr_squarem_iter(siris,
//                                                             sigma_beta(j),
//                                                             logodds(i),
//                                                             betahat.col(g),
//                                                             se.col(g),
//                                                             copy_alpha,
//                                                             copy_mu,
//                                                             copy_SiRiSr,
//                                                             tolerance,
//                                                             itermax,
//                                                             lnz_tol,iter,max_err);
//                      sigbvec(mat_ind)=sigma_beta(j);
//                      logoddsvec(mat_ind)=logodds(i);
//                      fgeneidmat(mat_ind)=fgeneid(g);
//                      nlzvec(mat_ind)=retvec;
//                      itervec(mat_ind)=iter;
//                      errvec(mat_ind)=max_err;
//                      npivec(mat_ind)=copy_alpha.mean();
//                      pvevec(mat_ind)=copy_SiRiSr.matrix().transpose()*(copy_alpha*copy_mu).matrix();
//                    }
//                  }
//                });
//   
//   
//   //  Rcpp::Rcout<<"mean lnZ is: "<<mean(nlzvec)<<std::endl;
//   return(Rcpp::DataFrame::create(_["logodds"]=logoddsvec,
//                                  _["sigb"]=sigbvec,
//                                  _["rel_err"]=errvec,
//                                  _["iterations"]=itervec,
//                                  _["alpha_mean"]=npivec,
//                                  _["pve"]=pvevec,
//                                  _["lnZ"]=nlzvec,
//                                  _["fgeneid"]=fgeneidmat));
// }

#else


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
  
  using namespace Rcpp;
  size_t sigb_size= sigma_beta.size();
  size_t logodds_size=logodds.size();
  size_t tot_size=sigb_size;
  if(tot_size!=logodds_size){
    Rcpp::stop("Length of sigma_beta must equal length of logodds");
  }
  
  Eigen::ArrayXd npivec(tot_size);
  Eigen::ArrayXd pvevec(tot_size);
  Eigen::ArrayXd nlzvec(tot_size);
  Eigen::ArrayXd sigbvec(tot_size);
  Eigen::ArrayXd lovec(tot_size);
  Eigen::ArrayXd errvec(tot_size);
  Eigen::ArrayXi itervec(tot_size);
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
    nlzvec[t]=retvec;
    itervec[t]=iter;
    errvec[t]=max_err;
    npivec[t]=copy_alpha.mean();
    pvevec[t]=copy_SiRiSr.matrix().transpose()*(copy_alpha*copy_mu).matrix();
  }
  return(Rcpp::DataFrame::create(_["logodds"]=logodds,
                                 _["sigb"]=sigma_beta,
                                 _["rel_err"]=errvec,
                                 _["iterations"]=itervec,
                                 _["alpha_mean"]=npivec,
                                 _["pve"]=pvevec,
                                 _["lnZ"]=nlzvec));
}

#endif



#if RCPP_PARALLEL_USE_TBB

using namespace tbb;


template<typename T> Rcpp::DataFrame grid_rss_varbvsr_naive(
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
  Eigen::ArrayXd errvec(tot_size);
  Eigen::ArrayXi itervec(tot_size);
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
                   //                   sigbvec(t)=sigma_beta(j);
                   Eigen::ArrayXd copy_alpha(talpha0);
                   Eigen::ArrayXd copy_mu(tmu0);
                   Eigen::ArrayXd copy_SiRiSr(tSiRiSr0);
                   //                   lovec(t)=logodds(i);
                   int iter=0;
                   double max_err=1;
                   double retvec=rss_varbvsr_naive_iter(SiRiS,
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
                   itervec[t]=iter;
                   errvec[t]=max_err;
                   npivec[t]=copy_alpha.mean();
                   pvevec[t]=copy_SiRiSr.matrix().transpose()*(copy_alpha*copy_mu).matrix();
                 }});
  
  //  Rcpp::Rcout<<"mean lnZ is: "<<mean(nlzvec)<<std::endl;
  return(Rcpp::DataFrame::create(_["logodds"]=logodds,
                                 _["sigb"]=sigma_beta,
                                 _["rel_err"]=errvec,
                                 _["iterations"]=itervec,
                                 _["alpha_mean"]=npivec,
                                 _["pve"]=pvevec,
                                 _["lnZ"]=nlzvec));
}




#else


template<typename T> Rcpp::DataFrame grid_rss_varbvsr_naive(
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
  Eigen::ArrayXd nlzvec(tot_size);
  Eigen::ArrayXd sigbvec(tot_size);
  Eigen::ArrayXd lovec(tot_size);
  Eigen::ArrayXd errvec(tot_size);
  Eigen::ArrayXi itervec(tot_size);
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
    
    
    double retvec =rss_varbvsr_naive_iter(SiRiS,
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
    itervec[t]=iter;
    errvec[t]=max_err;
    npivec[t]=copy_alpha.mean();
    pvevec[t]=copy_SiRiSr.matrix().transpose()*(copy_alpha*copy_mu).matrix();
  }
  return(Rcpp::DataFrame::create(_["logodds"]=logodds,
                                 _["sigb"]=sigma_beta,
                                 _["rel_err"]=errvec,
                                 _["iterations"]=itervec,
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
    const Rcpp::NumericMatrix SiRiS,
    const Rcpp::NumericVector sigma_beta,
    const Rcpp::NumericVector logodds,
    const Rcpp::NumericVector betahat,
    const Rcpp::NumericVector se,
    const Rcpp::NumericVector talpha0,
    const Rcpp::NumericVector tmu0,
    const Rcpp::NumericVector tSiRiSr0,
    double tolerance,
    int itermax,
    Rcpp::LogicalVector verbose,
    Rcpp::LogicalVector lnz_tol){
  //  std::cout<<"Starting grid_search_rss_varbvsr"<<std::endl;
  bool isVerbose = verbose(0);
  bool islnz_tol = lnz_tol(0);
  typedef typename RcppParallel::RMatrix<double> PRmat;
  typedef typename RcppParallel::RVector<double> PRvec;
  
  
 
  const PRmat _SiRiS(SiRiS);
  const PRvec _sigma_beta(sigma_beta);
  const PRvec _logodds(logodds);
  const PRvec _betahat(betahat);
  const PRvec _se(se);
  const PRvec _talpha0(talpha0);
  const PRvec _tmu0(tmu0);
  const PRvec _tSiRiSr0(tSiRiSr0);
  
  const c_Matrix_external __SiRiS(_SiRiS.begin(),SiRiS.nrow(),SiRiS.ncol());
  const c_arrayxd_external __sigma_beta(_sigma_beta.begin(),sigma_beta.size());
  const c_arrayxd_external __logodds(_logodds.begin(),logodds.size());
  
  const c_arrayxd_external __betahat(_betahat.begin(),betahat.size());
  const c_arrayxd_external __se(_se.begin(),se.size());
  const c_arrayxd_external __talpha0(_talpha0.begin(),talpha0.size());
  const c_arrayxd_external __tmu0(_tmu0.begin(),tmu0.size());
  const c_arrayxd_external __tSiRiSr0(_tSiRiSr0.begin(),tSiRiSr0.size());
  if(__sigma_beta.minCoeff()<=0){
    Rcpp::stop("sigma_beta must be strictly positive");
  }
  return grid_rss_varbvsr<c_Matrix_internal>(__SiRiS,__sigma_beta,__logodds,__betahat,
                                             __se,__talpha0,__tmu0,__tSiRiSr0,tolerance,
                                             itermax,isVerbose,islnz_tol);
} 



// 
// Rcpp::DataFrame grid_search_rss_varbvsr_multitrait(
//     const  Matrix_external R,
//     const arrayxd_external sigma_beta,
//     const arrayxd_external logodds,
//     const Matrix_external  betahat,
//     const Matrix_external  se,
//     const arrayxd_external talpha0,
//     const arrayxd_external tmu0,
//     const arrayxi_external fgeneid,
//     double tolerance,
//     int itermax,
//     Rcpp::LogicalVector verbose,
//     Rcpp::LogicalVector lnz_tol){
//   //  std::cout<<"Starting grid_search_rss_varbvsr"<<std::endl;
//   bool isVerbose = verbose(0);
//   bool islnz_tol = lnz_tol(0);
//   
//   if(sigma_beta.minCoeff()<=0){
//     Rcpp::stop("sigma_beta must be strictly positive");
//   }
//   
//   return grid_rss_varbvsr_multitrait(R,sigma_beta,logodds,betahat,
//                                              se,talpha0,tmu0,tolerance,
//                                              itermax,isVerbose,fgeneid,islnz_tol);
// } 





//[[Rcpp::export]]
Rcpp::DataFrame grid_search_rss_varbvsr_naive(
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
  
  return grid_rss_varbvsr_naive<c_Matrix_internal>(SiRiS,sigma_beta,logodds,betahat,
                                                   se,talpha0,tmu0,tSiRiSr0,tolerance,
                                                   itermax,isVerbose,islnz_tol);
} 






//[[Rcpp::export]]
Rcpp::DataFrame grid_search_rss_varbvsr_naive_sp(
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
  
  return grid_rss_varbvsr_naive<c_sparseMatrix_internal>(SiRiS,sigma_beta,logodds,betahat,
                                                   se,talpha0,tmu0,tSiRiSr0,tolerance,
                                                   itermax,isVerbose,islnz_tol);
} 







