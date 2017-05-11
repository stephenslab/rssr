#include <RcppEigen.h>
#include "rssr.h"
#include <cstdio>
#include <cmath>
#include <RcppParallel.h>

#include <progress.hpp>
// [[Rcpp::depends(RcppProgress)]]
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
                                             double &mtp, double &lnZ,const double lnZ0,const bool reverse)
{
  
  
  // alpha_r=alpha1-alpha0;
  // alpha_v=(alpha-alpha1)-alpha_r;
  // 
  // mu_r   = mu1-mu0;
  // mu_v   = mu-mu1-mu_r;
//  Rcpp::Rcerr<<"Do we need to backtrack?"<<std::endl;
  double sigma_beta_square=sigma_beta*sigma_beta;
  
  
  if((mtp<(-1)) && (lnZ < lnZ0)){
    size_t num_bt=0;
//    Rcpp::Rcerr<<"Starting backtracking"<<std::endl;
    while((lnZ<lnZ0) &&(num_bt < 10)){
      mtp = 0.5*(mtp-1);
      alpha = alpha0-2*mtp*(alpha1-alpha0)+(mtp*mtp)*(alpha-2*alpha1+alpha0);
      mu = mu0-2*mtp*(mu1-mu0)+(mtp*mtp)*(mu-2*mu1+mu0);
      SiRiSr = (SiRiS*(alpha*mu).matrix()).array();
      rss_varbvsr_iter(SiRiS,sigma_beta_square,s,logodds,betahat,sesquare,ssrat,alpha,mu,SiRiSr,reverse);
      // r=alpha*mu;
      lnZ=calculate_lnZ(betahat/sesquare,alpha*mu,SiRiSr,logodds,sesquare,alpha,mu,s,sigma_beta);          
      num_bt=num_bt+1;
    }
    // if(num_bt>=10){
    //   Rcpp::Rcerr<<"backtrack failed"<<std::endl;
    // }
  }
}



// 
// template <typename T> void squarem_backtrack(arrayxxd_internal alpha,
//                                              const c_arrayxxd_internal alpha0,
//                                              const c_arrayxxd_internal alpha1,
//                                              arrayxxd_internal mu,
//                                              const c_arrayxxd_internal mu0,
//                                              const c_arrayxxd_internal mu1,
//                                              arrayxxd_internal SiRiSr,
//                                              T SiRiS,
//                                              const c_arrayxd_internal sigma_beta,
//                                              const c_arrayxd_internal logodds,
//                                              const c_arrayxxd_internal s,
//                                              const c_arrayxd_internal betahat,
//                                              const c_arrayxd_internal sesquare,
//                                              const c_arrayxxd_internal ssrat,
//                                              arrayxd_internal mtp, arrayxd_internal lnZ,
//                                              const c_arrayxd_internal lnZ0,
//                                              const bool reverse)
// {
//   
//   size_t tot_size=sigma_beta.size();
//   
//   for(size_t c=0; c<tot_size; c++){
//     //    std::cout<<"c:"<<c<<std::endl;
//     //  lnZ(c)=calculate_lnZ(betahat/sesquare,alpha.col(c)*mu.col(c),SiRiSr.col(c),logodds(c),sesquare,alpha.col(c),mu.col(c),s.col(c),sigma_beta(c));
//     //      std::cout<<"lnZ(c):"<<lnZ(c)<<std::endl;
//     // if(!lnZ.allFinite()){
//     //   Rcpp::stop("lnZ(c) is NaN");
//     // }
//     //      lnZ.coeff(i)=  calculate_lnZ(q,alpha*mu,SiRiSr,logodds,sesquare,alpha,mu,s,sigma_beta);
//     double tmtp = mtp.coeff(c);
//     
//     // alpha_r=alpha1-alpha0;
//     // alpha_v=(alpha-alpha1)-alpha_r;
//     // 
//     // mu_r   = mu1-mu0;
//     // mu_v   = mu-mu1-mu_r;
//     
//     
//     if((tmtp<(-1)) && (lnZ.coeff(c) < lnZ0.coeff(c))){
//       // Rcpp::Rcout<<"Begin bt"<<std::endl;
//       size_t num_bt=0;
//       while((lnZ.coeff(c)<lnZ0.coeff(c)) &&(num_bt < 10)){
//         
//         double sigma_beta_square=sigma_beta(c)*sigma_beta(c);
//         double tlogodds=logodds(c);
//         tmtp = 0.5*(tmtp-1);
//         
//         alpha.col(c) = alpha0.col(c)-2*tmtp*(alpha1-alpha0).col(c)+(tmtp*tmtp)*(alpha-2*alpha1+alpha0).col(c);
//         mu.col(c) = mu0.col(c)-2*tmtp*(mu1-mu0).col(c)+(tmtp*tmtp)*(mu-2*mu1+mu0).col(c);
//         SiRiSr.col(c) = SiRiS*(alpha.col(c)*mu.col(c)).matrix();
//         rss_varbvsr_iter(SiRiS,sigma_beta_square,s.col(c),tlogodds,betahat,sesquare,ssrat.col(c),alpha.col(c),mu.col(c),SiRiSr.col(c),reverse);
//         lnZ(c)=calculate_lnZ(betahat/sesquare,alpha.col(c)*mu.col(c),SiRiSr.col(c),logodds(c),sesquare,alpha.col(c),mu.col(c),s.col(c),sigma_beta(c));
//         num_bt=num_bt+1;
//       }
//     }
//   }
// }


double calculate_mtp(const c_arrayxd_internal alpha,const c_arrayxd_internal alpha1,const c_arrayxd_internal alpha0, const c_arrayxd_internal mu,const c_arrayxd_internal mu1, const c_arrayxd_internal mu0){
  return(  -std::sqrt(((alpha1-alpha0).square()).sum()+((mu1-mu0).square()).sum())/std::sqrt(((alpha-2*alpha1+alpha0).square()).sum()+((mu-2*mu1+mu0).square()).sum()));
}

double calculate_mtp(const c_arrayxd_internal alpha_r,const c_arrayxd_internal alpha_v,const c_arrayxd_internal mu_r, const c_arrayxd_internal mu_v){
  
  return(  -std::sqrt(((alpha_r).square()).sum()+((mu_r).square()).sum())/std::sqrt(((alpha_v).square()).sum()+((mu_v).square()).sum()));
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
  // Rcpp::Rcout<<"mtp::"<<mtp<<std::endl;
  // Rcpp::Rcout<<"tmtp::"<<tmtp<<std::endl;
  // if((alpha*mu).hasNaN()){
  //   Rcpp::Rcerr<<" (pre-adjustment)"<<std::endl;
  //   Rcpp::stop("alpha*mu is not finite!");
  // }
  // 
  // if(alpha_r.hasNaN()){
  //   Rcpp::stop("alpha_r is not finite");
  // }
  // if(alpha_v.hasNaN()){
  //   Rcpp::stop("alpha_v is not finite");
  // }
  // if((alpha0).hasNaN()){
  //   Rcpp::stop("alpha0 is not finite");
  // }
  // if((2*mtp*(alpha_r)).hasNaN()){
  //   Rcpp::Rcerr<<mtp<<"_"<<2*mtp<<std::endl;
  //   Rcpp::Rcerr<<"alpha_r"<<alpha_r.square().sum()<<std::endl;
  //   Rcpp::Rcerr<<"alpha_v"<<alpha_v.square().sum()<<std::endl;
  //   Rcpp::Rcerr<<alpha_r.minCoeff()<<"_"<<2*mtp*(alpha_r.minCoeff())<<std::endl;
  //   Rcpp::Rcerr<<alpha_r.maxCoeff()<<"_"<<2*mtp*(alpha_r.maxCoeff())<<std::endl;
  //   Rcpp::stop("alpha update_1_2 is not finite");
  // }
  // if((mtpsq*(alpha_v)).hasNaN()){
  //   Rcpp::stop("alpha update_2 is not finite");
  // }
  
  

  if(mtp >=-1){
  }else{
   if(!std::isnan(mtp)){
    alpha=alpha0-2*mtp*(alpha_r)+(mtp*mtp)*(alpha_v);
    mu=mu0-2*mtp*(mu_r)+(mtp*mtp)*(mu_v);
    SiRiSr=(SiRiS*(alpha*mu).matrix()).array();
  }
  }
  
}






// template <typename T> void squarem_adjust(arrayxxd_internal alpha,
//                                           const c_arrayxxd_internal alpha0,
//                                           const c_arrayxxd_internal alpha1,
//                                           arrayxxd_internal mu,
//                                           const c_arrayxxd_internal mu0,
//                                           const c_arrayxxd_internal mu1,
//                                           arrayxxd_internal SiRiSr,
//                                           T SiRiS,
//                                           arrayxd_internal mtp){
//   
//   size_t tot_size=mtp.size();
//   
//   for(size_t c=0; c<tot_size; c++){
//     double tmtp= -sqrt(((alpha1-alpha0).col(c).square()).sum()+((mu1-mu0).col(c).square()).sum())/sqrt(((alpha-2*alpha1+alpha0).col(c).square()).sum()+((mu-2*mu1+mu0).col(c).square()).sum());
//     mtp(c)=tmtp;
//     if(tmtp >=-1){
//     }else{
//       alpha.col(c)=alpha0.col(c)-2*tmtp*(alpha1-alpha0).col(c)+(tmtp*tmtp)*(alpha-2*alpha1+alpha0).col(c);
//       mu.col(c)=mu0.col(c)-(mu1-mu0).col(c)*2*tmtp+(mu-2*mu1+mu0).col(c)*(tmtp*tmtp);
//       SiRiSr.col(c)=SiRiS*(alpha.col(c)*mu.col(c)).matrix();
//     }
//   }
// }


// void calc_max_err(const c_arrayxd_internal lnZ,const c_arrayxd_internal lnZ0,const c_arrayxxd_internal alpha, const c_arrayxxd_internal alpha0,const c_arrayxxd_internal mu,const c_arrayxxd_internal mu0,  double &max_err,const bool lnztol){
//   size_t tot_size=lnZ.size();
//   if(lnztol){
//     max_err=((lnZ-lnZ0).abs()/(lnZ.abs()+lnZ0.abs()+double_lim::epsilon())).maxCoeff();
//   }else{
//     max_err=0;
//     for(size_t c=0; c<tot_size; c++){
//       double tmax=find_maxerr(alpha.col(c),alpha0.col(c),alpha.col(c)*mu.col(c),alpha0.col(c)*mu0.col(c));
//       if(tmax>max_err){
//         max_err=tmax;
//       }
//     }
//   }
// }

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

// Eigen::ArrayXd initialize_value(const double val, const c_arrayxd_internal param){
//   size_t tot_size= param.size();
//   Eigen::ArrayXd ret(tot_size);
//   for(size_t i=0; i<tot_size;i++){
//     ret(i)=val;
//   }
//   return(ret);
// }

Eigen::ArrayXd initialize_array(const c_arrayxd_internal init, const double param){
  Eigen::ArrayXd ret=init;
  return(ret);
}

// Eigen::ArrayXXd initialize_array(const c_arrayxd_internal init, const c_arrayxd_internal param){
//   size_t p=init.size();
//   size_t tot_size= param.size();
//   Eigen::ArrayXXd ret(p,tot_size);
//   for(size_t i=0; i<tot_size;i++){
//     ret.col(i)=init;
//   }
//   return(ret);
// }

Eigen::ArrayXd initialize_s(const c_arrayxd_internal sesquare, const double sigma_beta_square){
  return((sesquare*(sigma_beta_square))/(sesquare+(sigma_beta_square)));
}

// Eigen::ArrayXXd initialize_s(const c_arrayxd_internal sesquare, const c_arrayxd_internal sigma_beta_square){
//   size_t p=sesquare.size();
//   size_t tot_size= sigma_beta_square.size();
//   Eigen::ArrayXXd ret(p,tot_size);
//   for(size_t i=0; i<tot_size;i++){
//     ret.col(i)=(sesquare*(sigma_beta_square(i)))/(sesquare+(sigma_beta_square(i)));
//   }
//   return(ret);
// }

Eigen::ArrayXd initialize_ssrat(const c_arrayxd_internal s, const double sigma_beta_square){
  return((s/sigma_beta_square).log());
}

// Eigen::ArrayXXd initialize_ssrat(const c_arrayxxd_internal s, const c_arrayxd_internal sigma_beta_square){
//   size_t p=s.rows();
//   size_t tot_size= sigma_beta_square.size();
//   Eigen::ArrayXXd ret(p,tot_size);
//   for(size_t i=0; i<tot_size;i++){
//     ret.col(i)=(s.col(i)/sigma_beta_square(i)).log();
//   }
//   return(ret);
// }

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
  //This function implements RSS with variational bayes and the SQUAREM algorithm.
  int iter=0;
  bool lnztol=lnz_tol[0];
  const size_t p = betahat.size();
  
  // Eigen::ArrayXd alpha=talpha0;
  // Eigen::ArrayXd mu=tmu0;
  // Eigen::ArrayXd SiRiSr=tSiRiSr0;
  double lnZ= 0;
  //  double lnZ=log(0);
  
  
  Eigen::ArrayXd alpha0 = alpha;
  Eigen::ArrayXd mu0 =mu;
  Eigen::ArrayXd  SiRiSr0 = SiRiSr;
  // Eigen::ArrayXd alpha0=alpha;
  // Eigen::ArrayXd mu0=mu;
  // Eigen::ArrayXd SiRiSr0=SiRiSr;
  
  
  ArrayXd alpha1=alpha;
  ArrayXd mu1=mu;
  ArrayXd SiRiSr1=SiRiSr;
  
  
  ArrayXd r=alpha*mu;
  
  
  double sigma_beta_square=sigma_beta*sigma_beta;
  Eigen::ArrayXd sesquare =se*se;
  Eigen::ArrayXd  q= betahat/sesquare;
  ArrayXd s=initialize_s(sesquare,sigma_beta_square);
  ArrayXd ssrat=initialize_ssrat(s,sigma_beta_square);
  // Eigen::ArrayXd  s= (sesquare*(sigma_beta*sigma_beta))/(sesquare+(sigma_beta*sigma_beta));
  // Eigen::ArrayXd ssrat((s/sigma_beta_square).log());
  
  
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
  

    

    Rcpp::Rcout<<"squarem adjust"<<std::endl;
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




// // 
// Rcpp::List wrap_squarem_backtrack(const arrayxd_external alpha,
//                                   const arrayxd_external alpha0,
//                                   const arrayxd_external alpha1,
//                                   const arrayxd_external mu,
//                                   const arrayxd_external mu0,
//                                   const arrayxd_external mu1,
//                                   const arrayxd_external SiRiSr,
//                                   const Matrix_external SiRiS,
//                                   const double sigma_beta,
//                                   const double logodds,
//                                   double mtp,
//                                   const arrayxd_external betahat,
//                                   const arrayxd_external se){
//   
//   
//   ArrayXd talpha=alpha;
//   ArrayXd tmu=mu;
//   ArrayXd tSiRiSr=SiRiSr;
//   double sigma_beta_square=sigma_beta*sigma_beta;
//   ArrayXd sesquare=se.square();
//   ArrayXd s=initialize_s(sesquare,sigma_beta_square);
//   ArrayXd ssrat=initialize_ssrat(s,sigma_beta_square);
// 
//   Eigen::ArrayXd  q= betahat/sesquare;
//   
//   double lnZ=  calculate_lnZ(q,talpha*tmu,SiRiSr,logodds,sesquare,talpha,tmu,s,sigma_beta);
//   
//   // if((mtp<(-1)) && (lnZ < lnZ0)){
//   //   size_t num_bt=0;
//   //   while((lnZ<lnZ0) &&(num_bt < 10)){
//   //     mtp = 0.5*(mtp-1);
//   //     alpha = alpha0-2*mtp*alpha_r+(mtp*mtp)*alpha_v;
//   //     mu = mu0-2*mtp*mu_r+(mtp*mtp)*mu_v;
//   //     SiRiSr = SiRiS*(alpha*mu).matrix();
//   //     rss_varbvsr_iter(SiRiS,sigma_beta_square,s,logodds,betahat,sesquare,ssrat,alpha,mu,SiRiSr,reverse);
//   //     r=alpha*mu;
//   //     lnZ=calculate_lnZ(q,r,SiRiSr,logodds,sesquare,alpha,mu,s,sigma_beta);          
//   //     num_bt=num_bt+1;
//   //   }
//   // }
//   
//   squarem_backtrack(talpha,alpha0,alpha1,
//                        tmu,mu0,mu1,
//                        tSiRiSr,SiRiS,
//                        sigma_beta,logodds,
//                        s,betahat,
//                        sesquare,ssrat,
//                        mtp,lnZ, lnZ0, reverse);
// 
//   
// }
                               

// 
// Rcpp::List wrap_squarem_adjust(const Matrix_external SiRiS,
//                                     const double sigma_beta,
//                                     const double logodds,
//                                     const arrayxd_external  betahat,
//                                     const arrayxd_external se,
//                                     arrayxd_external talpha,
//                                     arrayxd_external tmu,
//                                     arrayxd_external tSiRiSr,
//                                     const double tolerance,
//                                     const int itermax,
//                                     Rcpp::LogicalVector lnz_tol){


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
  double lnZ= initialize_value(0,logodds);
  //  double lnZ=log(0);
  
  
  Eigen::ArrayXd alpha0 = alpha;
  Eigen::ArrayXd mu0 =mu;
  Eigen::ArrayXd  SiRiSr0 = SiRiSr;
  // Eigen::ArrayXd alpha0=alpha;
  // Eigen::ArrayXd mu0=mu;
  // Eigen::ArrayXd SiRiSr0=SiRiSr;
  
  
  ArrayXd alpha1=alpha;
  ArrayXd mu1=mu;
  ArrayXd SiRiSr1=SiRiSr;
  
  
  ArrayXd r=alpha*mu;
  
  
  double sigma_beta_square=sigma_beta*sigma_beta;
  Eigen::ArrayXd sesquare =se*se;
  Eigen::ArrayXd  q= betahat/sesquare;
  ArrayXd s=initialize_s(sesquare,sigma_beta_square);
  ArrayXd ssrat=initialize_ssrat(s,sigma_beta_square);
  // Eigen::ArrayXd  s= (sesquare*(sigma_beta*sigma_beta))/(sesquare+(sigma_beta*sigma_beta));
  // Eigen::ArrayXd ssrat((s/sigma_beta_square).log());
  
  ArrayXd alpha_r=(alpha1-alpha0);
  ArrayXd alpha_v=(alpha-alpha1)-alpha_r;
  
  ArrayXd mu_r=(mu1-mu0);
  ArrayXd mu_v=(mu-mu1)-mu_r;
  if(r.hasNaN()){
    Rcpp::Rcerr<<"In iteration iter 0(0)"<<std::endl;
    Rcpp::stop("alpha*mu is not finite!");
  }
  
  lnZ=calculate_lnZ(q,r,SiRiSr,logodds,sesquare,alpha,mu,s,sigma_beta);
  
  double mtp=0;
  
  max_err=1;
  double lnZ0=lnZ;
  
  
  while(max_err>tolerance){
    lnZ0=lnZ;
    alpha0=alpha;
    mu0=mu; 
    
    bool reverse = iter%2!=0;
    // if((alpha*mu).hasNaN()){
    //   Rcpp::Rcerr<<"alpha*mu is not finite In iteration iter(0) : "<<iter<<std::endl;
    //   //      Rcpp::stop("alpha*mu is not finite!");
    // }//else{
    //   Rcpp::Rcerr<<"alpha*mu IS finite In iteration iter(0) : "<<iter<<std::endl;
    // }
    rss_varbvsr_iter(SiRiS,sigma_beta_square,s,logodds,betahat,sesquare,ssrat,alpha,mu,SiRiSr,reverse);
    
    alpha1=alpha;
    mu1=mu;
    SiRiSr1=SiRiSr;
    
    // if((alpha*mu).hasNaN()){
    //   Rcpp::Rcerr<<"alpha*mu is not finite In iteration iter(1) : "<<iter<<std::endl;
    //   //      Rcpp::stop("alpha*mu is not finite!");
    // }//else{
    //   Rcpp::Rcerr<<"alpha*mu IS finite In iteration iter(1) : "<<iter<<std::endl;
    // }
    rss_varbvsr_iter(SiRiS,sigma_beta_square,s,logodds,betahat,sesquare,ssrat,alpha,mu,SiRiSr,reverse);
    
    
    
    // if((alpha*mu).hasNaN()){
    //   Rcpp::Rcerr<<"alpha*mu is not finite In iteration iter (2.1): "<<iter<<std::endl;
    //   //      Rcpp::stop("alpha*mu is not finite!");
    // }//else{
    // //   Rcpp::Rcerr<<"alpha*mu IS finite In iteration iter(2.1) : "<<iter<<std::endl;
    // }
    
    alpha_r=(alpha1-alpha0);
    alpha_v=(alpha-alpha1)-alpha_r;
    
    mu_r=(mu1-mu0);
    mu_v=(mu-mu1)-mu_r;
    squarem_adjust<T>(alpha,alpha0,alpha1,alpha_r,alpha_v,mu,mu0,mu1,mu_r,mu_v,SiRiSr,SiRiS,mtp);
    // if((alpha*mu).hasNaN()){
    //   Rcpp::Rcerr<<"alpha*mu is not finite In iteration iter (2.2): "<<iter<<std::endl;
    //   //      Rcpp::stop("alpha*mu is not finite!");
    // }//else{
    //   Rcpp::Rcerr<<"alpha*mu IS finite In iteration iter(2.2) : "<<iter<<std::endl;
    // }
    
    rss_varbvsr_iter(SiRiS,sigma_beta_square,s,logodds,betahat,sesquare,ssrat,alpha,mu,SiRiSr,reverse);
    //    r=alpha*mu;
    // if((alpha*mu).hasNaN()){
    //   Rcpp::Rcerr<<"alpha*mu is not finite In iteration iter (3): "<<iter<<std::endl;
    //   //      Rcpp::stop("alpha*mu is not finite!");
    //  }//else{
    //   Rcpp::Rcerr<<"alpha*mu IS finite In iteration iter(3) : "<<iter<<std::endl;
    // }
    
    lnZ=  calculate_lnZ(q,alpha*mu,SiRiSr,logodds,sesquare,alpha,mu,s,sigma_beta);

    squarem_backtrack<T>(alpha,alpha0,alpha1,
                         mu,mu0,mu1,
                         SiRiSr,SiRiS,
                         sigma_beta,logodds,
                         s,betahat,
                         sesquare,ssrat,
                         mtp,lnZ, lnZ0, reverse);
    if((alpha*mu).hasNaN()){
      Rcpp::Rcerr<<"alpha*mu is not finite In iteration iter (4): "<<iter<<std::endl;
      //      Rcpp::stop("alpha*mu is not finite!");
    }//else{
    //   Rcpp::Rcerr<<"alpha*mu IS finite In iteration iter(4) : "<<iter<<std::endl;
    // }
    
    calc_max_err(lnZ,lnZ0, alpha, alpha0, mu, mu0,max_err,lnztol);
    
    if(iter>itermax){
      printf("Maximum iteration number reached: %+0.2d \n",(int)iter);

      printf("The log variational lower bound of the last step increased by(at most) %+0.2e\n",max_err);
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
  Rcpp::NumericVector nlzvec(tot_size);
  Rcpp::NumericVector sigbvec(tot_size);
  Rcpp::NumericVector lovec(tot_size);
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
  Rcpp::NumericVector nlzvec(tot_size);
  Rcpp::NumericVector sigbvec(tot_size);
  Rcpp::NumericVector lovec(tot_size);
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


// 
// 
// using namespace tbb;
// 
// 
// template<typename T> Rcpp::DataFrame grid_rss_varbvsr_array(
//     const T SiRiS,
//     const c_arrayxd_internal sigma_beta,
//     const c_arrayxd_internal logodds,
//     const c_arrayxd_internal betahat,
//     const c_arrayxd_internal se,
//     const c_arrayxd_internal talpha0,
//     const c_arrayxd_internal tmu0,
//     const c_arrayxd_internal tSiRiSr0,
//     double tolerance,
//     int itermax,
//     bool isVerbose,
//     bool islnz_tol){
//   
//   //  std::cout<<"Starting grid_rss_varbvsr (tbb)"<<std::endl;
//   
//   using namespace Rcpp;
//   size_t sigb_size= sigma_beta.size();
//   size_t logodds_size=logodds.size();
//   
//   size_t tot_size=sigb_size;
//   if(tot_size!=logodds_size){
//     Rcpp::stop("Length of sigma_beta must equal length of logodds");
//   }
//   
//   Eigen::ArrayXd npivec(tot_size);
//   Eigen::ArrayXd nlzvec(tot_size);
//   // Eigen::ArrayXd sigbvec(tot_size);
//   // Eigen::ArrayXd lovec(tot_size);
//   Eigen::ArrayXd pvevec(tot_size);
//   Eigen::ArrayXd errvec(tot_size);
//   Eigen::ArrayXi itervec(tot_size);
//   
//   
//   Rcpp::LogicalVector verbose(1);
//   verbose(0)=isVerbose;
//   Rcpp::LogicalVector lnz_tol(1);
//   lnz_tol(0)=islnz_tol;
//   
//   parallel_for(blocked_range<size_t>(0,tot_size),
//                [&](const blocked_range<size_t>& r)  {
//                  Eigen::ArrayXd sub_sigma_beta=sigma_beta.segment(r.begin(),r.end());
// //                 sigbvec.segment(r.begin(),r.end())=sub_sigma_beta;
//                  Eigen::ArrayXd sub_logodds=logodds.segment(r.begin(),r.end());
// //                 lovec.segment(r.begin(),r.end())=sub_logodds;
//                  Eigen::ArrayXXd copy_alpha = initialize_array(talpha0,sub_sigma_beta);
//                  Eigen::ArrayXXd copy_mu = initialize_array(tmu0,sub_sigma_beta);
//                  Eigen::ArrayXXd copy_SiRiSr = initialize_array(tSiRiSr0,sub_sigma_beta);
//                  int iter=0;
//                  double max_err=1;
//                  Eigen::ArrayXd retvec=rss_varbvsr_squarem_iter(SiRiS,
//                                                                 sub_sigma_beta,
//                                                                 sub_logodds,
//                                                                 betahat,
//                                                                 se,
//                                                                 copy_alpha,
//                                                                 copy_mu,
//                                                                 copy_SiRiSr,
//                                                                 tolerance,
//                                                                 itermax,
//                                                                 lnz_tol,
//                                                                 iter,max_err);
//                  // for(size_t t=r.begin(); t!=r.end(); t++){
//                  //   size_t i=t;
//                  //   size_t j=t;
//                  //   sigbvec(t)=sigma_beta(j);
//                  //   Eigen::ArrayXd copy_alpha(talpha0);
//                  //   Eigen::ArrayXd copy_mu(tmu0);
//                  //   Eigen::ArrayXd copy_SiRiSr(tSiRiSr0);
//                  //   lovec(t)=logodds(i);
//                  //   int iter=0;
//                  //   double max_err=1;
//                  //   double retvec=rss_varbvsr_squarem_iter(SiRiS,
//                  //                                          sigma_beta(j),
//                  //                                          logodds(i),
//                  //                                          );
// 
//                    for(int i=r.begin(); i!=r.end();i++){
//                      nlzvec(i)=retvec(i-r.begin());
//                      itervec(i)=iter;
//                      errvec(i)=max_err;
//                      npivec(i)=copy_alpha.col(i-r.begin()).mean();
//                      pvevec(i)=copy_SiRiSr.col(i-r.begin()).matrix().transpose()*(copy_alpha.col(i-r.begin())*copy_mu.col(i-r.begin())).matrix();
//                    }
// 
// //                   npivec.segment(r.begin(),r.end())=copy_alpha.rowwise().mean();
// //                   pvevec.segment(r.begin(),r.end())=(copy_SiRiSr.matrix().transpose()*(copy_alpha*copy_mu).matrix()).array();
//                  });
//   
//   //  Rcpp::Rcout<<"mean lnZ is: "<<mean(nlzvec)<<std::endl;
//   return(Rcpp::DataFrame::create(_["logodds"]=logodds,
//                                  _["sigb"]=sigma_beta,
//                                  _["rel_err"]=errvec,
//                                  _["iterations"]=itervec,
//                                  _["alpha_mean"]=npivec,
//                                  _["pve"]=pvevec,
//                                  _["lnZ"]=nlzvec));
// }
// 
// 
// //[[Rcpp::export]]
// Rcpp::DataFrame grid_search_rss_varbvsr_array(
//     const Rcpp::NumericMatrix SiRiS,
//     const Rcpp::NumericVector sigma_beta,
//     const Rcpp::NumericVector logodds,
//     const Rcpp::NumericVector  betahat,
//     const Rcpp::NumericVector  se,
//     const Rcpp::NumericVector talpha0,
//     const Rcpp::NumericVector tmu0,
//     const Rcpp::NumericVector tSiRiSr0,
//     double tolerance,
//     int itermax,
//     Rcpp::LogicalVector verbose,
//     Rcpp::LogicalVector lnz_tol){
//   //  std::cout<<"Starting grid_search_rss_varbvsr"<<std::endl;
//   bool isVerbose = verbose(0);
//   bool islnz_tol = lnz_tol(0);
// 
//   using namespace Eigen;
//   
//   MatrixXd tSiRiS = MatrixXd::Map(SiRiS.begin(),SiRiS.nrow(),SiRiS.ncol());
//   
//   ArrayXd tsigma_beta= ArrayXd::Map(sigma_beta.begin(),sigma_beta.size());
//   ArrayXd tlogodds= ArrayXd::Map(logodds.begin(),logodds.size());
//   
//   ArrayXd tbetahat= ArrayXd::Map(betahat.begin(),betahat.size());
//   ArrayXd tse= ArrayXd::Map(se.begin(),se.size());
//   
//   ArrayXd alpha0= ArrayXd::Map(talpha0.begin(),talpha0.size());
//   ArrayXd mu0= ArrayXd::Map(tmu0.begin(),tmu0.size());
//   ArrayXd SiRiSr0= ArrayXd::Map(tSiRiSr0.begin(),tSiRiSr0.size());
//   
//   
// 
//   
//   if(tsigma_beta.minCoeff()<=0){
//     Rcpp::stop("sigma_beta must be strictly positive");
//   }
//   
//   return grid_rss_varbvsr_array<c_Matrix_internal>(tSiRiS,tsigma_beta,tlogodds,tbetahat,
//                                              tse,alpha0,mu0,SiRiSr0,tolerance,
//                                              itermax,isVerbose,islnz_tol);
// } 




