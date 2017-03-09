#include <RcppEigen.h>
#include "kl.hpp"
#include "rss_varbvsr_iter_alt.hpp"
#include "sigmoid.hpp"
#include <math.h>








//[[Rcpp::export]]
Rcpp::List rss_varbvsr_naive_alt (const Eigen::SparseMatrix<double> R_matrix,
                              const double sigma_beta,
                              const double logodds,
                              const Eigen::ArrayXd betahat,
                              const Eigen::ArrayXd se,
                              const Eigen::ArrayXd &talpha0,
                              const Eigen::ArrayXd &tmu0,
                              const Eigen::ArrayXd &tSiRiSr0,
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
  double lnZ=log(0);
  Eigen::MatrixXd params0(p,2);
  Eigen::MatrixXd params(p,2);
  
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
    bool reverse = iter%2!=0;
    rss_varbvsr_iter_alt(R_matrix,sigma_beta,logodds,betahat,sesquare,alpha,mu,SiRiSr,reverse);
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



//' Run RSS with the variational bayes algorithm accelerated with SQUAREM
//' @template rssr
//' @param talpha0 a length p vector specifying the initial value of alpha
//' @param tmu0 a length p vector specifying the initial value of mu
//' @param SiRiSr0 a length p vector specifying the initial value of SiRiSr
//[[Rcpp::export]]
Rcpp::List rss_varbvsr_squarem_alt(const Eigen::SparseMatrix<double> &SiRiS,
                                   const double sigma_beta,
                                   const double logodds,
                                   const Eigen::ArrayXd betahat,
                                   const Eigen::ArrayXd se,
                                   const Eigen::ArrayXd &talpha0,
                                   const Eigen::ArrayXd &tmu0,
                                   const Eigen::ArrayXd &tSiRiSr0,
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
    rss_varbvsr_iter_alt(SiRiS,sigma_beta,logodds,betahat,se,alpha,mu,SiRiSr,reverse);
    assert(alpha.min()>=0);
    alpha1=alpha;
    mu1=mu;
    SiRiSr1=SiRiSr;
    
    rss_varbvsr_iter_alt(SiRiS,sigma_beta,logodds,betahat,se,alpha,mu,SiRiSr,reverse);
    assert(alpha.min()>=0);
    
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
    assert(alpha.min()>=0);
    lnZ=  calculate_lnZ(q,alpha*mu,SiRiSr,logodds,sesquare,alpha,mu,s,sigma_beta);
    if((mtp<(-1)) && (lnZ < lnZ0)){
      size_t num_bt=0;
      while((lnZ<lnZ0) &&(num_bt < 10)){
        mtp = 0.5*(mtp-1);
        alpha = alpha0-2*mtp*alpha_r+(mtp*mtp)*alpha_v;
        assert(alpha.min()>=0);
        mu = mu0-2*mtp*mu_r+(mtp*mtp)*mu_v;
        SiRiSr = SiRiS*(alpha*mu).matrix();
        rss_varbvsr_iter_alt(SiRiS,sigma_beta,logodds,betahat,se,alpha,mu,SiRiSr,reverse);
        assert(alpha.min()>=0);
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













// 
// void rss_varbvsr_iter_alt(const Eigen::SparseMatrix<double> SiRiS,
//                           const double sigma_beta,
//                           const double logodds,
//                           const Eigen::ArrayXd betahat,
//                           const Eigen::ArrayXd se,
//                           Eigen::ArrayXd &alpha,
//                           Eigen::ArrayXd &mu,
//                           Eigen::ArrayXd &SiRiSr,
//                           bool reverse){
//   
//   
//   // Get the number of SNPs (p) and coordinate ascent updates (m).
//   const size_t p = betahat.size();
//   
//   Eigen::ArrayXd alpha0=alpha;
//   Eigen::ArrayXd mu0=mu;
//   Eigen::SparseMatrix<double> SiRiS0=SiRiS;
//   for(size_t i=0; i<p; i++){
//     SiRiS0.coeffRef(i,i)=0;
//   }
//   
//   Eigen::ArrayXd se_square = se * se;
//   double  sigma_beta_square = sigma_beta * sigma_beta;
//   
//   // Compute the variational estimate of the posterior variance.
//   Eigen::ArrayXd sigma_square = (se_square * sigma_beta_square) / (se_square + sigma_beta_square);
//   // Update the variational estimate of the posterior mean.
//   Eigen::ArrayXd r0 = alpha0 * mu0;
//   mu = sigma_square * (betahat / se_square - (SiRiS0*(r0.matrix())).array());
//   
//   alpha = (logodds + 0.5 * ((sigma_square/(sigma_beta_square)).log() + (mu0 * mu0) / sigma_square)).unaryExpr(std::ptr_fun(sigmoid));
//   
//   SiRiSr  =(SiRiS*(alpha*mu).matrix()).array();
//   
// }
















// USAGE: run a single coordinate ascent of variational update to fit RSS-BVSR

// void rss_varbvsr_update_alt (const double betahat,
//                          const double se,
//                          const double sigma_beta,
//                          const Eigen::ArrayXd &SiRiS_snp,
//                          Eigen::ArrayXd &SiRiSr,
//                          const double SiRiSr_snp,
//                          const double logodds,
//                          double &alpha,
//                          double &mu) {
//   
//  
// }



void rss_varbvsr_iter_alt(const Eigen::SparseMatrix<double> R_matrix,
                              const double sigma_beta,
                              const double logodds,
                              const Eigen::ArrayXd betahat,
                              const Eigen::ArrayXd sesquare,
                              Eigen::ArrayXd &alpha,
                              Eigen::ArrayXd &mu,
                              Eigen::ArrayXd &SiRiSr,
                              bool reverse){
  
  
  // Get the number of SNPs (p) and coordinate ascent updates (m).
  const size_t p = betahat.size();
  
  // Initialize outputs.
  
  // Store a single column of matrix inv(S)*R*inv(S).
  Eigen::ArrayXd alpha0=alpha;
  Eigen::ArrayXd mu0=mu;
//  Eigen::ArrayXd SiRiSr0=SiRiSr;
  double sigma_beta_square = sigma_beta * sigma_beta;
  Eigen::SparseMatrix<double> t_R_matrix=R_matrix;
  for(size_t t=0;t<p;t++){
    t_R_matrix.coeffRef(t,t)=0;
  }
  Eigen::ArrayXd  SiRiS_snp(p);
  Eigen::VectorXd  SiRiS_snp_v(p);
  // Run coordinate ascent updates.
  // Repeat for each coordinate ascent update.
  size_t i=0;
  for (size_t j = 0; j < p; j++) {
    if(reverse){
      i=p-1-j;
    }else{
      i=j;
    }
    SiRiS_snp_v=(tSiRiS.col(i));
    SiRiS_snp=SiRiS_snp_v.array();
    // Copy the kth column of matrix inv(S)*R*inv(S).
    // copySparseColumn(SiRiS_snp, k, SiRiS.elems, Ir, Jc, p);
    
    // Copy the kth element of vector inv(S)*R*inv(S)*r.
    double SiRiSr_snp = SiRiSr0(i);
    
    // Perform the mean-field variational update.
    
   
    
    // Compute the variational estimate of the posterior variance.
    double sigma_square = (sesquare(i) * sigma_beta_square) / (sesquare(i) + sigma_beta_square);
    
    // Update the variational estimate of the posterior mean.
    double r = alpha(i) * mu(i);
double tsum=0;
for(size_t t=0; t<)
    mu(i) = sigma_square * (betahat(i) / sesquare(i) - SiRiSr_snp);
    
    // Update the variational estimate of the posterior inclusion probability.
    double SSR = mu(i) * mu(i) / sigma_square;
    alpha(i) = sigmoid(logodds + 0.5 * (log(sigma_square/(sigma_beta_square)) + SSR));
    
    // Update SiRiSr = inv(S)*R*inv(S)*r
    double r_new = alpha(i) * mu(i);
    SiRiSr0+=(SiRiS_snp*(r_new-r));
   // rss_varbvsr_update_alt(betahat(i), se(i), sigma_beta, 
                           // SiRiS_snp, SiRiSr, SiRiSr_snp, 
                           // logodds, alpha(i), mu(i));
  }
  SiRiSr=(SiRiS*(alpha*mu).matrix()).array();
}








//[[Rcpp::export]]
Rcpp::List wrap_rss_varbvsr_iter_alt(const Eigen::SparseMatrix<double> R_matrix,
                                     const double sigma_beta,
                                     const double logodds,
                                     const Eigen::ArrayXd betahat,
                                     const Eigen::ArrayXd se,
                                     Eigen::ArrayXd &alpha,
                                     Eigen::ArrayXd &mu,
                                     Eigen::ArrayXd &SiRiSr,
                                     bool reverse){
  
  Eigen::ArrayXd sesquare=se*se;
  rss_varbvsr_iter_alt(R_matrix,
                       sigma_beta,
                       logodds,
                       betahat,
                       sesquare,
                       alpha,
                       mu,
                       SiRiSr,
                       reverse);
  return Rcpp::List::create(Rcpp::Named("alpha1")=alpha,
                            Rcpp::Named("mu1")=mu,
                            Rcpp::Named("SiRiSr")=SiRiSr);
}










