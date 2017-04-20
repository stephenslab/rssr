#include <RcppEigen.h>
#include "rssr.h"


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
  
  return grid_rss_varbvsr(SiRiS,sigma_beta,logodds,betahat,
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
  
  return grid_rss_varbvsr(SiRiS,sigma_beta,logodds,betahat,
                          se,talpha0,tmu0,tSiRiSr0,tolerance,
                          itermax,isVerbose,islnz_tol);
} 





//[[Rcpp::export]]
Rcpp::DataFrame grid_search_rss_varbvsr_array(const Matrix_external SiRiS,
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
  
  size_t sigb_size= sigma_beta.size();
  size_t logodds_size=logodds.size();
  size_t tot_size=sigb_size*logodds_size;
  
  Eigen::ArrayXd nlzvec(tot_size);
  
  Eigen::ArrayXd sigbvec(tot_size);
  Eigen::ArrayXd lovec(tot_size);
  for(int t=0; t<tot_size;t++){
    size_t i=t%logodds_size;
    size_t j=t/logodds_size;
    sigbvec(t)=sigma_beta(j);
    lovec(t)=logodds(i);
  }
  
  return(rss_varbvsr_squarem_array(SiRiS,sigbvec,lovec,betahat,se,talpha0,tmu0,tSiRiSr0,tolerance,itermax,verbose(0),lnz_tol(0)));
}



// 
// 
// Rcpp::DataFrame grid_search_rss_varbvsr_array(
//     const  Matrix_external SiRiS,
//     const arrayxd_external sigma_beta,
//     const arrayxd_external logodds,
//     const arrayxd_external  betahat,
//     const arrayxd_external  se,
//     const arrayxd_external talpha0,
//     const arrayxd_external tmu0,
//     const arrayxd_external tSiRiSr0,
//     double tolerance,
//     int itermax,
//     Rcpp::LogicalVector verbose,
//     Rcpp::LogicalVector lnz_tol){
// //  std::cout<<"Starting grid_search_rss_varbvsr"<<std::endl;
//   bool isVerbose = verbose(0);
//   bool islnz_tol = lnz_tol(0);
//   
//   return grid_rss_varbvsr_array(SiRiS,sigma_beta,logodds,betahat,
//                           se,talpha0,tmu0,tSiRiSr0,tolerance,
//                           itermax,isVerbose,islnz_tol);
// } 
// 


