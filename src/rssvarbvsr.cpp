#include <RcppEigen.h>
#include "rssr.h"

#include <math.h>

using namespace Rcpp;

 
void rss_varbvsr_update(const double betahat,
                        const double se_square,
                        const double sigma_beta_square,
                        const c_arrayxd_internal SiRiS_snp,
                        const double sigma_square,
                        arrayxd_internal SiRiSr,
                        const double SiRiSr_snp,
                        const double ssrat,
                        const double logodds,
                        double &alpha,
                        double &mu) {

  
  const size_t p=SiRiS_snp.size();

  
  double r = alpha * mu;
  
  mu = sigma_square * (betahat / se_square + r/se_square - SiRiSr_snp);
  alpha = sigmoid(logodds + 0.5 * (ssrat + mu*mu/sigma_square));
  
  double r_new = alpha * mu-r;
  SiRiSr+=SiRiS_snp*r_new;
  
}






void rss_varbvsr_iter(const c_Matrix_internal SiRiS,
                      const double sigma_beta_square,
                      const c_arrayxd_internal sigma_square,
                      const double logodds,
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
    rss_varbvsr_update(betahat.coeff(i),
                       se_square.coeff(i),
                       sigma_beta_square, 
                       SiRiS.col(i),
                       sigma_square.coeff(i),
                       SiRiSr,
                       SiRiSr.coeff(i),
                       ssrat.coeff(i),
                       logodds,
                       alpha.coeffRef(i),
                       mu.coeffRef(i));
  }
}




void rss_varbvsr_iter(const c_sparseMatrix_internal SiRiS,
                      const double sigma_beta_square,
                      const c_arrayxd_internal sigma_square,
                      const double logodds,
                      const c_arrayxd_internal betahat,
                      const c_arrayxd_internal se_square,
                      const c_arrayxd_internal ssrat,
                      arrayxd_internal alpha,
                      arrayxd_internal mu,
                      arrayxd_internal SiRiSr,
                      bool reverse){
  
  
  // Get the number of SNPs (p) and coordinate ascent updates (m).
  const size_t p = betahat.size();
  
  // Initialize outputs.
  
  // Store a single column of matrix inv(S)*R*inv(S).
  
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
    SiRiS_snp_v=(SiRiS.col(i));
    SiRiS_snp=SiRiS_snp_v.array();
    // Copy the kth column of matrix inv(S)*R*inv(S).
    // copySparseColumn(SiRiS_snp, k, SiRiS.elems, Ir, Jc, p);
    
    // Copy the kth element of vector inv(S)*R*inv(S)*r.
    double SiRiSr_snp = SiRiSr.coeff(i);
    
    
    
    rss_varbvsr_update(betahat.coeff(i),
                       se_square.coeff(i),
                       sigma_beta_square, 
                       SiRiS_snp,
                       sigma_square.coeff(i),
                       SiRiSr,
                       SiRiSr.coeff(i),
                       ssrat.coeff(i),
                       logodds,
                       alpha.coeffRef(i),
                       mu.coeffRef(i));
  }
}


// 
// void rss_varbvsr_iter(const c_Matrix_internal SiRiS,
//                        const c_arrayxd_internal sigma_beta_square,
//                        const c_arrayxxd_internal sigma_square,
//                        const c_arrayxd_internal logodds,
//                        const c_arrayxd_internal betahat,
//                        const c_arrayxd_internal se_square,
//                        const c_arrayxxd_internal ssrat,
//                        arrayxxd_internal alpha,
//                        arrayxxd_internal mu,
//                        arrayxxd_internal SiRiSr,
//                        bool reverse){
//   using namespace Rcpp;
//   size_t tot_size = sigma_beta_square.size();
//   // Get the number of SNPs (p) and coordinate ascent updates (m).
//   const size_t p = betahat.size();
//   
//   // Initialize outputs.
//   
//   // Store a single column of matrix inv(S)*R*inv(S).
//   
//   // Eigen::ArrayXd  SiRiS_snp(p);
//   // Eigen::VectorXd  SiRiS_snp_v(p);
//   
//   // Run coordinate ascent updates.
//   // Repeat for each coordinate ascent update.
//   size_t i=0;
//   
//   // Eigen::ArrayXd se_square=se*se;
//   // RowArray sigma_beta_square = (sigma_beta*sigma_beta).transpose();
//   // RowArray logodds=tlogodds.transpose();
//   //  Eigen::ArrayXd sigma_beta_square=sigma_beta*sigma_beta;
//   RowArray r(tot_size);
//   RowArray r_new(tot_size);
// 
// 
// 
//   
// 
//   
//   // mu = sigma_square * (betahat / se_square + r/se_square - SiRiSr_snp);
//   // alpha = sigmoid(logodds + 0.5 * (ssrat + mu*mu/sigma_square));
//   
//   // double r_new = alpha * mu-r;
//   // SiRiSr+=SiRiS_snp*r_new;
//   
//   // Update SiRiSr = inv(S)*R*inv(S)*r
//   
//   //  RowArray new_alpha(tot_size);
//   for (size_t j = 0; j < p; j++) {
//     if(reverse){
//       i=p-1-j;
//     }else{
//       i=j;
//     }
//     r = alpha.row(i) * mu.row(i);
//     mu.row(i) = sigma_square.row(i) * (betahat.coeff(i) / se_square.coeff(i) + r / se_square.coeff(i) - SiRiSr.row(i));
//     
//     // Update the variational estimate of the posterior inclusion probability.
// 
//     //    new_alpha = 1/(1+(-(logodds + ((sigma_square.row(i)/(sigma_beta_square)).log() + SSR)*0.5)).exp());
//     // if(j<2){
//     //   std::cout<<j<<": SSR_size"<<SSR.size()<<std::endl;
//     //   std::cout<<" (sigma_square.row(i)/sigma_beta_square): "<<(sigma_square.row(i)/sigma_beta_square).size()<<std::endl;
//     //   std::cout<<" (sigma_square.row(i)/sigma_beta_square).log(): "<<(sigma_square.row(i)/sigma_beta_square).log().size()<<std::endl;
//     //   std::cout<<" logodds.size(): "<<logodds.size();
//     //   std::cout<<"sigma_square.row(i).size(): "<<sigma_square.row(i).size();
//     //   std::cout<<" sigma_beta_square.size(): "<<sigma_beta_square.size();
//     //   std::cout<<" tot_size:"<<tot_size<<std::endl;
//     //   std::cout<<j<<": "<<new_alpha<<" :"<<new_alpha.size()<<std::endl;
//     // }
//     alpha.row(i) =1/(1+(-(logodds + (ssrat.row(i) + (mu.row(i) * mu.row(i)) / sigma_square.row(i))*0.5)).exp());
//     
//     // Update SiRiSr = inv(S)*R*inv(S)*r
//     r_new = alpha.row(i) * mu.row(i)-r;
//     for(size_t c=0; c<r_new.size(); c++){
//       SiRiSr.col(c)+=(SiRiS.col(i).array())*(r_new.coeff(c));
//     }
//   }
//   
//   
// }
// 





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
  double sigma_beta_square=sigma_beta*sigma_beta;
  Eigen::ArrayXd sesquare=se.square();
  Eigen::ArrayXd  s= (sesquare*(sigma_beta*sigma_beta))/(sesquare+(sigma_beta*sigma_beta));
  Eigen::ArrayXd ssrat((s/sigma_beta_square).log());
  rss_varbvsr_iter(SiRiS,sigma_beta_square,s,logodds,betahat,sesquare,ssrat,talpha,tmu,tSiRiSr,reverse);
  return Rcpp::List::create(Rcpp::Named("alpha1")=talpha,
                            Rcpp::Named("mu1")=tmu,
                            Rcpp::Named("SiRiSr")=tSiRiSr);
}

//[[Rcpp::export]]
Eigen::ArrayXd wrap_rss_varbvsr_iter_squarem(const arrayxd_external alpha_mu,
                                         const Matrix_external SiRiS,
                                         const double sigma_beta,
                                         const double logodds,
                                         const arrayxd_external betahat,
                                         const arrayxd_external se){
  
  
  size_t p=betahat.size();
  // const arrayxd_external alpha,
  // const arrayxd_external mu,
  // const arrayxd_external SiRiSr,
  using namespace Eigen;  
  
  Eigen::ArrayXd talpha=alpha_mu.head(p);
  Eigen::ArrayXd tmu=alpha_mu.tail(p);
  Eigen::ArrayXd tSiRiSr=SiRiS*(talpha*tmu).matrix();
  Eigen::ArrayXd ret_alpha_mu(alpha_mu.size());
  double sigma_beta_square=sigma_beta*sigma_beta;
  Eigen::ArrayXd sesquare=se.square();
  Eigen::ArrayXd  s= (sesquare*(sigma_beta*sigma_beta))/(sesquare+(sigma_beta*sigma_beta));
  Eigen::ArrayXd ssrat((s/sigma_beta_square).log());
  rss_varbvsr_iter(SiRiS,sigma_beta_square,s,logodds,betahat,sesquare,ssrat,talpha,tmu,tSiRiSr,true);
  rss_varbvsr_iter(SiRiS,sigma_beta_square,s,logodds,betahat,sesquare,ssrat,talpha,tmu,tSiRiSr,false);
  ret_alpha_mu<<talpha,tmu;
  return(ret_alpha_mu);
}






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
  double sigma_beta_square=sigma_beta*sigma_beta;
  Eigen::ArrayXd sesquare=se.square();
  Eigen::ArrayXd  s= (sesquare*(sigma_beta*sigma_beta))/(sesquare+(sigma_beta*sigma_beta));
  Eigen::ArrayXd ssrat((s/sigma_beta_square).log());
  rss_varbvsr_iter(SiRiS,sigma_beta_square,s,logodds,betahat,sesquare,ssrat,talpha,tmu,tSiRiSr,reverse);
  return Rcpp::List::create(Rcpp::Named("alpha1")=talpha,
                            Rcpp::Named("mu1")=tmu,
                            Rcpp::Named("SiRiSr")=tSiRiSr);
  
}











