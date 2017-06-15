#include <RcppEigen.h>
#include "rssr.h"

#include <math.h>
#include <cstdio>
#include <cmath>
#include <RcppParallel.h>

#include <progress.hpp>
// [[Rcpp::depends(RcppProgress)]]
//[[Rcpp::depends(RcppParallel)]]



using namespace Rcpp;
using namespace tbb;


//Try with intgamma=0?

double calculate_lnZ(const c_vectorxd_internal q,
                     const c_vectorxd_internal r,
                     const c_vectorxd_internal SiRiSr,
                     const c_vectorxd_internal sesquare,
                     const c_vectorxd_internal alpha,
                     const c_vectorxd_internal mu,
                     const c_vectorxd_internal s,
                     double sigb){
  
  
  double lnz0 = q.dot(r)-0.5*r.dot(SiRiSr);
  // if(!std::isfinite(lnz0)){
  //   Rcpp::stop("lnZ0 is not finite");
  // }
  double lnz1 = lnz0-0.5*(1/sesquare.array()).matrix().dot(betavar(alpha.array(),mu.array(),s.array()).matrix());
  // if(!std::isfinite(lnz0)){
  //   Rcpp::stop("lnZ1 is not finite");
  // }
  double lnz2 = lnz1+intklbeta_rssbvsr(alpha.array(),mu.array(),s.array(),sigb*sigb);
  // if(!std::isfinite(lnz0)){
  //   Rcpp::stop("lnZ2 is not finite");
  // }
  return(lnz2);
}


void rss_varbvsr_update(const double betahat,
                        const double se_square,
                        const double sigma_beta_square,
                        const c_arrayxd_internal SiRiS_snp,
                        const double sigma_square,
                        arrayxd_internal SiRiSr,
                        const double SiRiSr_snp,
                        const double ssrat,
                        double &alpha,
                        double &mu) {
  
  
  const size_t p=SiRiS_snp.size();
  
  
  double mu0 = mu;
  
  mu = sigma_square * (betahat / se_square + mu0/se_square - SiRiSr_snp);
  // alpha = sigmoid(logodds + 0.5 * (ssrat + mu*mu/sigma_square));
  
  double mu_new = mu-mu0;
  SiRiSr+=SiRiS_snp*mu_new;
  
}




void rss_varbvsr_iter(const c_Matrix_internal SiRiS,
                      const double sigma_beta_square,
                      const c_arrayxd_internal sigma_square,
                      const c_arrayxd_internal betahat,
                      const c_arrayxd_internal se_square,
                      const c_arrayxd_internal ssrat,
                      arrayxd_internal alpha,
                      arrayxd_internal mu,
                      arrayxd_internal SiRiSr,
                      bool reverse){
  
  // mkl_set_num_threads_local(1);
  size_t p=betahat.size();
  // Get the number of SNPs (p) and coordinate ascent updates (m).
  
  // Initialize outputs.
  
  // Store a single column of matrix inv(S)*R*inv(S).
  
  
  // Run coordinate ascent updates.
  // Repeat for each coordinate ascent update.
  
  size_t i=0;
  for (size_t j = 0; j < p; j++) {
    if(reverse){
      i=p-1-j;
    }else{
      i=j;
    }
    
    // Copy the kth column of matrix inv(S)*R*inv(S).
    // copySparseColumn(SiRiS_snp, k, SiRiS.elems, Ir, Jc, p);
    
    // Copy the kth element of vector inv(S)*R*inv(S)*r.
    
    
    // Perform the mean-field variational update.
    rss_varbvsr_update(betahat.coeffRef(i),
                       se_square.coeffRef(i),
                       sigma_beta_square, 
                       SiRiS.col(i),
                       sigma_square.coeffRef(i),
                       SiRiSr,
                       SiRiSr.coeffRef(i),
                       ssrat.coeffRef(i),
                       alpha.coeffRef(i),
                       mu.coeffRef(i));
  }
  
}



double calculate_mtp(const c_arrayxd_internal mu_r, const c_arrayxd_internal mu_v){
  
  return(  -std::sqrt(((mu_r).square()).sum())/std::sqrt(((mu_v).square()).sum()));
}

void squarem_adjust(arrayxd_internal alpha,
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
                    const c_Matrix_internal SiRiS,
                    double &mtp){
  
  
  
  mtp=calculate_mtp(mu_r,mu_v);
  
  
  
  
  if(mtp > (-1)){
    mtp=-1;
  }else{
    if(!std::isnan(mtp)){
      // alpha=alpha0-2*mtp*(alpha_r)+(mtp*mtp)*(alpha_v);
      mu=mu0-2*mtp*(mu_r)+(mtp*mtp)*(mu_v);
      SiRiSr=(SiRiS*(alpha*mu).matrix()).array();
    }
  }
}



template <typename T> void squarem_backtrack(arrayxd_internal &alpha,
                                             const c_arrayxd_internal alpha0,
                                             const c_arrayxd_internal alpha1,
                                             arrayxd_internal &mu,
                                             const c_arrayxd_internal mu0,
                                             const c_arrayxd_internal mu1,
                                             arrayxd_internal &SiRiSr,
                                             T SiRiS,
                                             const double sigma_beta,
                                             const c_arrayxd_internal s,
                                             const c_arrayxd_internal betahat,
                                             const c_arrayxd_internal sesquare,
                                             const c_arrayxd_internal ssrat,
                                             double &mtp, double &lnZ, double &lnZ0,double &lnZ00,bool& doDie)
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
     // alpha = alpha0-2*mtp*(alpha1-alpha0)+(mtp*mtp)*(alpha-2*alpha1+alpha0);
      mu = mu0-2*mtp*(mu1-mu0)+(mtp*mtp)*(mu-2*mu1+mu0);
      SiRiSr = (SiRiS*(alpha*mu).matrix()).array();
      rss_varbvsr_iter(SiRiS,sigma_beta_square,s,betahat,sesquare,ssrat,alpha,mu,SiRiSr,false);
      rss_varbvsr_iter(SiRiS,sigma_beta_square,s,betahat,sesquare,ssrat,alpha,mu,SiRiSr,true);
      // r=alpha*mu;
      lnZ=calculate_lnZ(betahat/sesquare,alpha*mu,SiRiSr,sesquare,alpha,mu,s,sigma_beta);          

      
      num_bt=num_bt+1;
      if(num_bt==max_backtrack){

        mu=mu0;
        SiRiSr = (SiRiS*(alpha*mu).matrix()).array();

        lnZ0=calculate_lnZ(betahat/sesquare,alpha*mu,SiRiSr,sesquare,alpha,mu,s,sigma_beta);
        rss_varbvsr_iter(SiRiS,sigma_beta_square,s,betahat,sesquare,ssrat,alpha,mu,SiRiSr,false);
        rss_varbvsr_iter(SiRiS,sigma_beta_square,s,betahat,sesquare,ssrat,alpha,mu,SiRiSr,true);
        lnZ=calculate_lnZ(betahat/sesquare,alpha*mu,SiRiSr,sesquare,alpha,mu,s,sigma_beta);
        if(lnZ<lnZ0){
          doDie=true;
        }
        mtp=-1;
      // Rcpp::Rcerr<<"num_bt:"<<num_bt<<" step_size:"<<mtp<<"lnZ:"<<lnZ<<" lnZ0:"<<lnZ0<<std::endl;
    }
  }
  }
  
}





template<typename T> double rss_varbvsr_squarem_iter(const T SiRiS,
                                                     const double sigma_beta,
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
  
  
  alpha.setOnes();
  SiRiSr=SiRiS*(alpha*mu).matrix();
  
  
  // Eigen::ArrayXd alpha=talpha0;
  // Eigen::ArrayXd mu=tmu0;
  // Eigen::ArrayXd SiRiSr=tSiRiSr0;
  double lnZ= initialize_value(0,1.0);
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
  
  lnZ=calculate_lnZ(q,r,SiRiSr,sesquare,alpha,mu,s,sigma_beta);
  
  
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
    
    rss_varbvsr_iter(SiRiS,sigma_beta_square,s,betahat,sesquare,ssrat,alpha,mu,SiRiSr,reverse);
    
    
    alpha1=alpha;
    mu1=mu;
    SiRiSr1=SiRiSr;
    
    
    rss_varbvsr_iter(SiRiS,sigma_beta_square,s,betahat,sesquare,ssrat,alpha,mu,SiRiSr,reverse);
    
    
    
    
    
    
    alpha_r=(alpha1-alpha0);
    alpha_v=(alpha-alpha1)-alpha_r;
    
    mu_r=(mu1-mu0);
    mu_v=(mu-mu1)-mu_r;
    squarem_adjust(alpha,alpha0,alpha1,alpha_r,alpha_v,mu,mu0,mu1,mu_r,mu_v,SiRiSr,SiRiS,mtp);
    
    
    rss_varbvsr_iter(SiRiS,sigma_beta_square,s,betahat,sesquare,ssrat,alpha,mu,SiRiSr,reverse);
    
    //    r=alpha*mu;
    // if((alpha*mu).hasNaN()){
    //   Rcpp::Rcerr<<"alpha*mu is not finite In iteration iter (3): "<<iter<<std::endl;
    //   //      Rcpp::stop("alpha*mu is not finite!");
    //  }//else{
    //   Rcpp::Rcerr<<"alpha*mu IS finite In iteration iter(3) : "<<iter<<std::endl;
    // }
    
    lnZ=  calculate_lnZ(q,alpha*mu,SiRiSr,sesquare,alpha,mu,s,sigma_beta);
    bool doDie=false;
    squarem_backtrack<T>(alpha,alpha0,alpha1,
                         mu,mu0,mu1,
                         SiRiSr,SiRiS,
                         sigma_beta,
                         s,betahat,
                         sesquare,ssrat,
                         mtp,lnZ, lnZ0,lnZ00, doDie);

    if((alpha*mu).hasNaN()){
      Rcpp::Rcerr<<"alpha*mu is not finite In iteration iter (4): "<<iter<<std::endl;
      //      Rcpp::stop("alpha*mu is not finite!");
    }//else{
    //   Rcpp::Rcerr<<"alpha*mu IS finite In iteration iter(4) : "<<iter<<std::endl;
    // }
    
    calc_max_err(lnZ,lnZ0, alpha, alpha0, mu, mu0,max_err,lnztol);
    if(doDie){
      iter=itermax+1;
    }
    
    if(iter>itermax){
//      printf("Maximum iteration number reached: %+0.2d \n",(int)iter);

//      printf("The log variational lower bound of the last step increased by(at most) %+0.2e\n",max_err);
      break;
    }
    iter=iter+1;
  }
  
  return lnZ;
}





template<typename T> double rss_varbvsr_naive_iter(const T SiRiS,
                                                   const double sigma_beta,
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
  double lnZ=0;
  alpha.setOnes();
  
  
  Eigen::ArrayXd talpha=alpha;
  Eigen::ArrayXd tmu=mu;
  Eigen::ArrayXd tSiRiSr=SiRiSr;
  
  double sigma_beta_square=sigma_beta*sigma_beta;
  Eigen::ArrayXd sesquare =se*se;
  Eigen::ArrayXd  q= betahat/sesquare;
  Eigen::ArrayXd  s= (sesquare*(sigma_beta*sigma_beta))/(sesquare+(sigma_beta*sigma_beta));
  Eigen::ArrayXd ssrat((s/sigma_beta_square).log());
  
  lnZ=calculate_lnZ(q,alpha*mu,SiRiSr,sesquare,alpha,mu,s,sigma_beta);
  
  
  max_err=1;
  double lnZ0=lnZ;
  double rel_l0=0;
  double rel_li=0;
  
  
  while(max_err>tolerance){
    talpha=alpha;
    tmu=mu;
    lnZ0=lnZ;
    bool reverse = iter%2!=0;
    rss_varbvsr_iter(SiRiS,sigma_beta_square,s,betahat,sesquare,ssrat,alpha,mu,SiRiSr,reverse);
    lnZ=  calculate_lnZ(q,alpha*mu,SiRiSr,sesquare,alpha,mu,s,sigma_beta);
    rel_li=rel_err(lnZ,lnZ0);
    if(lnztol){
      max_err=rel_li;
    }else{
      max_err=find_maxerr(alpha,talpha,alpha*mu,talpha*tmu);
    }
    
    iter=iter+1;
  }
  lnZ=  calculate_lnZ(q,alpha*mu,SiRiSr,sesquare,alpha,mu,s,sigma_beta);
  return(lnZ);
}






template<typename T> Rcpp::DataFrame grid_rss_varbvsr(
    const T SiRiS,
    const c_arrayxd_internal sigma_beta,
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
  
  size_t tot_size=sigb_size;

  
  Eigen::ArrayXd npivec(tot_size);

  Eigen::ArrayXd errvec(tot_size);
  Eigen::ArrayXi itervec(tot_size);
  Eigen::ArrayXd pvevec(tot_size);
//  pvevec.setZeros();
  Rcpp::NumericVector nlzvec(tot_size);
  
  Rcpp::LogicalVector verbose(1);
  verbose(0)=isVerbose;
  Rcpp::LogicalVector lnz_tol(1);
  lnz_tol(0)=islnz_tol;

  // static affinity_partitioner ap;

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
  return(Rcpp::DataFrame::create(_["sigb"]=sigma_beta,
                                 _["rel_err"]=errvec,
                                 _["iterations"]=itervec,
                                 _["alpha_mean"]=npivec,
                                 _["pve"]=pvevec,
                                 _["lnZ"]=nlzvec));
}





template<typename T> Rcpp::DataFrame grid_rss_varbvsr_naive(
    const T SiRiS,
    const c_arrayxd_internal sigma_beta,
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

  
  size_t tot_size=sigb_size;

  
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
  return(Rcpp::DataFrame::create(_["sigb"]=sigma_beta,
                                 _["rel_err"]=errvec,
                                 _["iterations"]=itervec,
                                 _["alpha_mean"]=npivec,
                                 _["pve"]=pvevec,
                                 _["lnZ"]=nlzvec));
}




//[[Rcpp::export]]
Rcpp::DataFrame grid_search_rss_varbvsr_norm(
    const  Matrix_external SiRiS,
    const arrayxd_external sigma_beta,
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
  
  return grid_rss_varbvsr<c_Matrix_internal>(SiRiS,sigma_beta,betahat,
                                             se,talpha0,tmu0,tSiRiSr0,tolerance,
                                             itermax,isVerbose,islnz_tol);
} 


//[[Rcpp::export]]
Rcpp::DataFrame grid_search_rss_varbvsr_naive_norm(
    const  Matrix_external SiRiS,
    const arrayxd_external sigma_beta,
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
  
  return grid_rss_varbvsr_naive<c_Matrix_internal>(SiRiS,sigma_beta,betahat,
                                                   se,talpha0,tmu0,tSiRiSr0,tolerance,
                                                   itermax,isVerbose,islnz_tol);
} 


//' Run RSS with the variational bayes algorithm accelerated with SQUAREM
//' @template rssr
//' @param talpha0 a length p vector specifying the initial value of alpha
//' @param tmu0 a length p vector specifying the initial value of mu
//' @param SiRiSr0 a length p vector specifying the initial value of SiRiSr
//[[Rcpp::export]]
Rcpp::List rss_varbvsr_squarem_norm(const Matrix_external SiRiS,
                               const double sigma_beta,
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
  
  double lnZ=rss_varbvsr_squarem_iter<c_Matrix_internal>(SiRiS,sigma_beta,betahat,se,alpha,mu,SiRiSr,tolerance,itermax,lnz_tol,iter,max_err);
  
  
  return Rcpp::List::create(Rcpp::Named("alpha")=alpha,
                            Rcpp::Named("mu")=mu,
                            Rcpp::Named("SiRiSr")=SiRiSr,
                            Rcpp::Named("max_err")=max_err,
                            Rcpp::Named("lnZ")=lnZ,
                            Rcpp::Named("iter")=iter);
}




//[[Rcpp::export]]
Rcpp::List rss_varbvsr_naive_norm (const Matrix_external SiRiS,
                              const double sigma_beta,
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
  
  double lnZ=rss_varbvsr_naive_iter<c_Matrix_internal>(SiRiS,sigma_beta,betahat,se,alpha,mu,SiRiSr,tolerance,itermax,lnz_tol,iter,max_err);
  
  return Rcpp::List::create(Rcpp::Named("alpha")=alpha,
                            Rcpp::Named("mu")=mu,
                            Rcpp::Named("SiRiSr")=SiRiSr,
                            Rcpp::Named("max_err")=max_err,
                            Rcpp::Named("lnZ")=lnZ,
                            Rcpp::Named("iter")=iter);
  
}
