#include <Rcpp.h>
#include <RcppEigen.h>
#include "rssr.h"
#include <math.h>
//[[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>

void rss_varbvsr_iter_alt(const c_Matrix_internal SiRiS,
                          const double sigma_beta,
                          const double logodds,
                          const c_arrayxd_internal betahat,
                          const c_arrayxd_internal se,
                          arrayxd_internal alpha,
                          arrayxd_internal mu,
                          arrayxd_internal SiRiSr,
                          bool reverse){
  
  
  // Get the number of SNPs (p) and coordinate ascent updates (m).
  const size_t p = betahat.size();
  
  // Initialize outputs.
  
  // Store a single column of matrix inv(S)*R*inv(S).
  
  // Eigen::ArrayXd  SiRiS_snp(p);
  // Eigen::VectorXd  SiRiS_snp_v(p);
  
  // Run coordinate ascent updates.
  // Repeat for each coordinate ascent update.
  Eigen::ArrayXd se_square=se*se;
  double sigma_beta_square = sigma_beta * sigma_beta;
  Eigen::ArrayXd sigma_square = (se_square * sigma_beta_square) / (se_square + sigma_beta_square);
  size_t i=0;
  for (size_t j = 0; j < p; j++) {
    if(reverse){
      i=p-1-j;
    }else{
      i=j;
    }
    // SiRiS_snp_v=(SiRiS.col(i));
    // SiRiS_snp=SiRiS_snp_v.array();
    // Copy the kth column of matrix inv(S)*R*inv(S).
    // copySparseColumn(SiRiS_snp, k, SiRiS.elems, Ir, Jc, p);
    
    // Copy the kth element of vector inv(S)*R*inv(S)*r.
    double SiRiSr_snp = SiRiSr.coeff(i);
    
    
  
    
    // Compute the variational estimate of the posterior variance.
    
    
    // Update the variational estimate of the posterior mean.
    double r = alpha.coeff(i) * mu.coeff(i);
    mu.coeffRef(i) = sigma_square.coeff(i) * (betahat.coeff(i) / se_square.coeff(i) + r / se_square.coeff(i) - SiRiSr_snp);
    
    // Update the variational estimate of the posterior inclusion probability.
    double SSR = mu.coeff(i) * mu.coeff(i) / sigma_square.coeff(i);
    alpha.coeffRef(i) = sigmoid(logodds + 0.5 * (log(sigma_square.coeff(i)/(sigma_beta_square)) + SSR));
    
    // Update SiRiSr = inv(S)*R*inv(S)*r
    double r_new = alpha.coeff(i) * mu.coeff(i);
    SiRiSr+=(SiRiS.col(i).array()*(r_new-r));
    
    // Perform the mean-field variational update.
    // rss_varbvsr_update(betahat(i), se(i), sigma_beta, 
    //                    SiRiS.col(i).array(), SiRiSr, SiRiSr_snp, 
    //                    logodds, alpha(i), mu(i));
  }
}

//[[Rcpp::export]]
Rcpp::List wrap_rss_varbvsr_iter_alt(const Matrix_external SiRiS,
                                      const double sigma_beta,
                                      const double logodds,
                                      const arrayxd_external betahat,
                                      const arrayxd_external se,
                                      const arrayxd_external alpha,
                                      const arrayxd_external mu,
                                      const arrayxd_external SiRiSr,
                                      Rcpp::LogicalVector reverse){

  Eigen::ArrayXd talpha=alpha;
  Eigen::ArrayXd tmu=mu;
  Eigen::ArrayXd tSiRiSr=SiRiSr;
  
  
  rss_varbvsr_iter_alt(SiRiS,
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
// 



//
//





double rss_varbvsr_squarem_iter_alt(const c_Matrix_internal SiRiS,
                                    const double sigma_beta,
                                    const double logodds,
                                    const c_arrayxd_internal betahat,
                                    const c_arrayxd_internal  se,
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
    rss_varbvsr_iter_alt(SiRiS,sigma_beta,logodds,betahat,se,alpha,mu,SiRiSr,reverse);
    alpha1=alpha;
    mu1=mu;
    SiRiSr1=SiRiSr;
    
    rss_varbvsr_iter_alt(SiRiS,sigma_beta,logodds,betahat,se,alpha,mu,SiRiSr,reverse);
    
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
    
    rss_varbvsr_iter_alt(SiRiS,sigma_beta,logodds,betahat,se,alpha,mu,SiRiSr,reverse);
    lnZ=  calculate_lnZ(q,alpha*mu,SiRiSr,logodds,sesquare,alpha,mu,s,sigma_beta);
    if((mtp<(-1)) && (lnZ < lnZ0)){
      size_t num_bt=0;
      while((lnZ<lnZ0) &&(num_bt < 10)){
        mtp = 0.5*(mtp-1);
        alpha = alpha0-2*mtp*alpha_r+(mtp*mtp)*alpha_v;
        mu = mu0-2*mtp*mu_r+(mtp*mtp)*mu_v;
        SiRiSr = SiRiS*(alpha*mu).matrix();
        rss_varbvsr_iter_alt(SiRiS,sigma_beta,logodds,betahat,se,alpha,mu,SiRiSr,reverse);
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







#if RCPP_PARALLEL_USE_TBB

using namespace tbb;


Rcpp::DataFrame grid_rss_varbvsr_alt(    const c_Matrix_internal SiRiS,
                                         c_arrayxd_internal sigma_beta,
                                         c_arrayxd_internal logodds,
                                         const c_arrayxd_internal betahat,
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

  Rcpp::LogicalVector verbose(1);
  verbose(0)=isVerbose;
  Rcpp::LogicalVector lnz_tol(1);
  lnz_tol(0)=islnz_tol;
  std::cout<<"Starting grid_rss_varbvsr_parallel (tbb)"<<std::endl;
  
  parallel_for(blocked_range<size_t>(0,tot_size),
               [&](const blocked_range<size_t>& r){
                 for(size_t t=r.begin(); t!=r.end(); t++){
                   size_t i=t%logodds_size;
                   size_t j=t/logodds_size;
                   sigbvec(t)=sigma_beta(j);
                   lovec(t)=logodds(i);
                   nlzvec(t)=rss_varbvsr_squarem_iter_alt(SiRiS,
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
  
  
  return(Rcpp::DataFrame::create(_["logodds"]=lovec,
                                 _["sigb"]=sigbvec,
                                 _["lnZ"]=nlzvec));
}




#else


Rcpp::DataFrame grid_rss_varbvsr_alt(const c_Matrix_internal SiRiS,
                                     const c_arrayxd_internal sigma_beta,
                                     const c_arrayxd_internal logodds,
                                     const c_arrayxd_internal betahat,
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

  std::cout<<"Starting grid_rss_varbvsr_alt (serial)"<<std::endl;
  Rcpp::LogicalVector verbose(1);
  verbose(0)=isVerbose;
  Rcpp::LogicalVector lnz_tol(1);
  lnz_tol(0)=islnz_tol;
  
  
  for(size_t t=0; t<tot_size; t++){
    size_t i=t%logodds_size;
    size_t j=t/logodds_size;
    sigbvec(t)=sigma_beta(j);
    lovec(t)=logodds(i);
    nlzvec(t)=rss_varbvsr_squarem_iter_alt(SiRiS,
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


//[[Rcpp::export]]
Rcpp::DataFrame grid_search_rss_varbvsr_alt(
    const Matrix_external SiRiS,
    const arrayxd_external sigma_beta,
    const arrayxd_external logodds,
    const arrayxd_external betahat,
    const arrayxd_external se,
    const arrayxd_external talpha0,
    const arrayxd_external tmu0,
    const arrayxd_external tSiRiSr0,
    double tolerance,
    int itermax,
    Rcpp::LogicalVector verbose,
    Rcpp::LogicalVector lnz_tol){
  
  bool isVerbose = verbose(0);
  bool islnz_tol = lnz_tol(0);
  return grid_rss_varbvsr_alt(SiRiS,sigma_beta,logodds,betahat,
                          se,talpha0,tmu0,tSiRiSr0,tolerance,
                          itermax,isVerbose,islnz_tol);
}  



