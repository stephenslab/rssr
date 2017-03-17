#include <Rcpp.h>
#include <RcppEigen.h>
#include "rssr.h"
#include <math.h>






// void rss_varbvsr_iter_alt_grid(const c_sparseMatrix_internal SiRiS,
//                                const c_arrayxd_internal sigma_beta,
//                                const c_arrayxd_internal tlogodds,
//                                const c_arrayxd_internal betahat,
//                                const c_arrayxd_internal se,
//                                arrayxxd_internal alpha,
//                                arrayxxd_internal mu,
//                                arrayxxd_internal SiRiSr,
//                                bool reverse){
//   
//   using namespace Rcpp;
//     size_t tot_size = sigma_beta.size();
//     // Get the number of SNPs (p) and coordinate ascent updates (m).
//     const size_t p = betahat.size();
//   for(size_t i=0; i<tot_size; i++){
//     rss_varbvsr_iter(SiRiS,sigma_beta(i),tlogodds(i),betahat,se,alpha.col(i),mu.col(i),SiRiSr.col(i),reverse);
//     
//   }
// }
  
  

void rss_varbvsr_iter_alt_grid(const c_sparseMatrix_internal SiRiS,
                               const c_arrayxd_internal sigma_beta,
                               const c_arrayxd_internal tlogodds,
                               const c_arrayxd_internal betahat,
                               const c_arrayxd_internal se,
                               const arrayxxd_internal sigma_square,
                               arrayxxd_internal alpha,
                               arrayxxd_internal mu,
                               arrayxxd_internal SiRiSr,
                               bool reverse){
  using namespace Rcpp;
  size_t tot_size = sigma_beta.size();
  // Get the number of SNPs (p) and coordinate ascent updates (m).
  const size_t p = betahat.size();

  // Initialize outputs.

  // Store a single column of matrix inv(S)*R*inv(S).

  // Eigen::ArrayXd  SiRiS_snp(p);
  // Eigen::VectorXd  SiRiS_snp_v(p);

  // Run coordinate ascent updates.
  // Repeat for each coordinate ascent update.
  size_t i=0;

  Eigen::ArrayXd se_square=se*se;
  RowArray sigma_beta_square = (sigma_beta*sigma_beta).transpose();
  RowArray logodds=tlogodds.transpose();
  //  Eigen::ArrayXd sigma_beta_square=sigma_beta*sigma_beta;
  RowArray r(tot_size);
  RowArray SSR(tot_size);
  RowArray r_new(tot_size);
  Eigen::VectorXd SiRiS_snp(p);
  // Update SiRiSr = inv(S)*R*inv(S)*r

  RowArray new_alpha(tot_size);
  for (size_t j = 0; j < p; j++) {
    if(reverse){
      i=p-1-j;
    }else{
      i=j;
    }
    r = alpha.row(i) * mu.row(i);
    mu.row(i) = sigma_square.row(i) * (betahat.coeff(i) / se_square.coeff(i) + r / se_square.coeff(i) - SiRiSr.row(i));

    // Update the variational estimate of the posterior inclusion probability.
    SSR = (mu.row(i) * mu.row(i)) / sigma_square.row(i);
    new_alpha = 1/(1+(-(logodds + ((sigma_square.row(i)/(sigma_beta_square)).log() + SSR)*0.5)).exp());
    // if(j<2){
    //   std::cout<<j<<": SSR_size"<<SSR.size()<<std::endl;
    //   std::cout<<" (sigma_square.row(i)/sigma_beta_square): "<<(sigma_square.row(i)/sigma_beta_square).size()<<std::endl;
    //   std::cout<<" (sigma_square.row(i)/sigma_beta_square).log(): "<<(sigma_square.row(i)/sigma_beta_square).log().size()<<std::endl;
    //   std::cout<<" logodds.size(): "<<logodds.size();
    //   std::cout<<"sigma_square.row(i).size(): "<<sigma_square.row(i).size();
    //   std::cout<<" sigma_beta_square.size(): "<<sigma_beta_square.size();
    //   std::cout<<" tot_size:"<<tot_size<<std::endl;
    //   std::cout<<j<<": "<<new_alpha<<" :"<<new_alpha.size()<<std::endl;
    // }
    alpha.row(i) =1/(1+(-(logodds + ((sigma_square.row(i)/(sigma_beta_square)).log() + SSR)*0.5)).exp());

    // Update SiRiSr = inv(S)*R*inv(S)*r
    r_new = alpha.row(i) * mu.row(i);
    SiRiS_snp=SiRiS.col(i);
    for(size_t c=0; c<r_new.size(); c++){

      SiRiSr.col(c)+=(SiRiS_snp.array())*(r_new.coeff(c)-r.coeff(c));
    }
  }


}


//[[Rcpp::export]]
Rcpp::List wrap_rss_varbvsr_iter_grid(const sparseMatrix_external SiRiS,
                                      const arrayxd_external sigma_beta,
                                      const arrayxd_external logodds,
                                      const arrayxd_external betahat,
                                      const arrayxd_external se,
                                      const arrayxd_external alpha,
                                      const arrayxd_external mu,
                                      const arrayxd_external SiRiSr,
                                      Rcpp::LogicalVector reverse){
  
  size_t p = betahat.size();
  size_t tot_size=sigma_beta.size();
  Eigen::ArrayXd sesquare=se*se;
  Eigen::ArrayXXd talpha(p,tot_size);
  talpha.colwise()=alpha;
  
  Eigen::ArrayXXd tmu(p,tot_size);
  tmu.colwise()=mu;
  Eigen::ArrayXXd tSiRiSr(p,tot_size);
  tSiRiSr.colwise()=SiRiSr;
  
  Eigen::ArrayXXd s(p,tot_size);
  for(size_t i=0; i<tot_size; i++){
    s.col(i)=(sesquare*(sigma_beta(i)*sigma_beta(i)))/(sesquare+(sigma_beta(i)*sigma_beta(i)));
  }
  
  rss_varbvsr_iter_alt_grid(SiRiS,
                            sigma_beta,
                            logodds,
                            betahat,
                            se,
                            s,
                            talpha,
                            tmu,
                            tSiRiSr,
                            reverse(0));
  return Rcpp::List::create(Rcpp::Named("alpha1")=talpha,
                            Rcpp::Named("mu1")=tmu,
                            Rcpp::Named("SiRiSr")=tSiRiSr);
}



Rcpp::DataFrame rss_varbvsr_alt_grid(const c_sparseMatrix_internal SiRiS,
                                     const c_arrayxd_internal sigma_beta,
                                     const c_arrayxd_internal logodds,
                                     const c_arrayxd_internal betahat,
                                     const c_arrayxd_internal se,
                                     const c_arrayxd_internal talpha0,
                                     const c_arrayxd_internal tmu0,
                                     const c_arrayxd_internal tSiRiSr0,
                                     double tolerance,
                                     int itermax,
                                     bool verbose,
                                     bool lnztol){
  
  
  
  //  bool lnztol=lnz_tol[0];
  using namespace Rcpp;
  using namespace Eigen;
  size_t sigb_size= sigma_beta.size();
  size_t logodds_size=logodds.size();
  size_t tot_size=sigb_size;
  if(tot_size!=logodds_size){
    Rcpp::stop("tot_size!=logodds_size");
  }
  size_t p=betahat.size();
  
  ArrayXXd alpha(p,tot_size);
  alpha.colwise()=talpha0;
  ArrayXXd mu(p,tot_size);
  mu.colwise()=tmu0;
  
  ArrayXXd SiRiSr(p,tot_size);
  SiRiSr.colwise()=tSiRiSr0;
  
  RowArray lnZ(tot_size);
  
  
  Eigen::ArrayXXd alpha0=alpha;
  Eigen::ArrayXXd mu0=mu;
//  Eigen::ArrayXXd SiRiSr0=SiRiSr;
  
  Eigen::ArrayXXd alpha1=alpha;
  Eigen::ArrayXXd mu1=mu;
//  Eigen::ArrayXXd SiRiSr1=SiRiSr;
  
  // Eigen::ArrayXXd alpha3=alpha;
  // Eigen::ArrayXXd mu3=mu;
//  Eigen::ArrayXXd SiRiSr3=SiRiSr;
  
  Eigen::ArrayXXd alpha_r(p,tot_size);
  Eigen::ArrayXXd alpha_v(p,tot_size);
  
  Eigen::ArrayXXd mu_v(p,tot_size);
  Eigen::ArrayXXd mu_r(p,tot_size);
  
  
  Eigen::ArrayXd sesquare =se*se;
  Eigen::ArrayXd  q= betahat/sesquare;
  Eigen::ArrayXXd s(p,tot_size);
//  Rcout<<"Initializing sigma_square and initial lnZ values"<<std::endl;
  for(size_t i=0; i<tot_size; i++){
    s.col(i)=(sesquare*(sigma_beta(i)*sigma_beta(i)))/(sesquare+(sigma_beta(i)*sigma_beta(i)));
    lnZ(i)=calculate_lnZ(q,alpha.col(i)*mu.col(i),SiRiSr.col(i),logodds(i),sesquare,alpha.col(i),mu.col(i),s.col(i),sigma_beta(i));
  }
  
  RowArray mtp(tot_size);
  int iter=0;
  //  iter.setZero();
  RowArray max_err(tot_size);
  max_err.setOnes();
  RowArray lnZ0=lnZ;

  
  RowArray rel_li(tot_size);
  rel_li.setZero();
//  Rcout<<"Starting main loop"<<std::endl;
  while(max_err.maxCoeff()>tolerance){
    //    Rcout<<iter<<": worst error is "<<max_err.maxCoeff()<<" tolerance is "<<tolerance<<std::endl;
    lnZ0=lnZ;
    alpha0=alpha;
    mu0=mu;
    
    
    
    bool reverse = iter%2!=0;
//    Rcout<<"First update"<<std::endl;
    rss_varbvsr_iter_alt_grid(SiRiS,sigma_beta,logodds,betahat,se,s,alpha,mu,SiRiSr,reverse);
    alpha1=alpha;
    mu1=mu;
//    SiRiSr1=SiRiSr;
//    Rcout<<"Second update"<<std::endl;
    rss_varbvsr_iter_alt_grid(SiRiS,sigma_beta,logodds,betahat,se,s,alpha,mu,SiRiSr,reverse);
    
    alpha_r=alpha1-alpha0;
    alpha_v=(alpha-alpha1)-alpha_r;
    
    mu_r   = mu1-mu0;
    mu_v   = mu-mu1-mu_r;
  
//    mtp= -sqrt((alpha_r.square()).colwise().sum()+(mu_r.square()).colwise().sum())/sqrt((alpha_v.square()).colwise().sum()+(mu_v.square()).colwise().sum());

//    Rcout<<"Calculating mtp"<<std::endl;
    for(size_t c=0; c<tot_size; c++){
      double tmtp= -sqrt((alpha_r.col(c).square()).sum()+(mu_r.col(c).square()).sum())/sqrt((alpha_v.col(c).square()).sum()+(mu_v.col(c).square()).sum());
      mtp(c)=tmtp;
      if(tmtp >=-1){
        
      }else{
        // Rcpp::Rcout<<"mtp <1: "<<tmtp<<std::endl;
        // Rcpp::Rcout<<"mean(alpha.col(c)): "<<alpha.col(c).mean()<<std::endl;
        // Rcpp::Rcout<<"mean(mu.col(c)): "<<mu.col(c).mean()<<std::endl;
        // Rcpp::Rcout<<"mean(SiRiSr.col(c)): "<<SiRiSr.col(c).mean()<<std::endl;
        alpha.col(c)=alpha0.col(c)-2*tmtp*alpha_r.col(c)+(tmtp*tmtp)*alpha_v.col(c);
        mu.col(c)=mu0.col(c)-mu_r.col(c)*2*tmtp+mu_v.col(c)*(tmtp*tmtp);
        SiRiSr.col(c)=SiRiS*(alpha.col(c)*mu.col(c)).matrix();
        // Rcpp::Rcout<<"mean(alpha.col(c)): "<<alpha.col(c).mean()<<std::endl;
        // Rcpp::Rcout<<"mean(mu.col(c)): "<<mu.col(c).mean()<<std::endl;
        // Rcpp::Rcout<<"mean(SiRiSr.col(c)): "<<SiRiSr.col(c).mean()<<std::endl;
      }
    }
//    Rcout<<"Third update"<<std::endl;
    
    rss_varbvsr_iter_alt_grid(SiRiS,sigma_beta,logodds,betahat,se,s,alpha,mu,SiRiSr,reverse);
//    Rcout<<"Calculating lower bound"<<std::endl;
    for(size_t c=0; c<tot_size; c++){
  //    std::cout<<"c:"<<c<<std::endl;
      lnZ(c)=calculate_lnZ(q,alpha.col(c)*mu.col(c),SiRiSr.col(c),logodds(c),sesquare,alpha.col(c),mu.col(c),s.col(c),sigma_beta(c));
//      std::cout<<"lnZ(c):"<<lnZ(c)<<std::endl;
      if(!lnZ.allFinite()){
        Rcpp::stop("lnZ(c) is NaN");
      }
      //      lnZ.coeff(i)=  calculate_lnZ(q,alpha*mu,SiRiSr,logodds,sesquare,alpha,mu,s,sigma_beta);
      double tmtp = mtp.coeff(c);
      
      
      if((tmtp<(-1)) && (lnZ.coeff(c) < lnZ0.coeff(c))){
        // Rcpp::Rcout<<"Begin bt"<<std::endl;
        size_t num_bt=0;
        while((lnZ.coeff(c)<lnZ0.coeff(c)) &&(num_bt < 10)){
          
          tmtp = 0.5*(tmtp-1);
          
          alpha.col(c) = alpha0.col(c)-2*tmtp*alpha_r.col(c)+(tmtp*tmtp)*alpha_v.col(c);
          mu.col(c) = mu0.col(c)-2*tmtp*mu_r.col(c)+(tmtp*tmtp)*mu_v.col(c);
          SiRiSr.col(c) = SiRiS*(alpha.col(c)*mu.col(c)).matrix();
          rss_varbvsr_iter(SiRiS,sigma_beta(c),logodds(c),betahat,se,alpha.col(c),mu.col(c),SiRiSr.col(c),reverse);
          lnZ(c)=calculate_lnZ(q,alpha.col(c)*mu.col(c),SiRiSr.col(c),logodds(c),sesquare,alpha.col(c),mu.col(c),s.col(c),sigma_beta(c));
          num_bt=num_bt+1;
        }
      }
    }
    //    return fabs(p0-p1)/(fabs(p0)+fabs(p1)+double_lim::epsilon());
    
    //    rel_l0=rel_err(lnZ00,lnZ);
    rel_li=(lnZ-lnZ0).abs()/(lnZ.abs()+lnZ0.abs()+double_lim::epsilon());
    if(lnztol){
      max_err=rel_li;
    }else{
      for(size_t c=0; c<tot_size; c++){
        max_err(c)=find_maxerr(alpha.col(c),alpha0.col(c),alpha.col(c)*mu.col(c),alpha0.col(c)*mu0.col(c));
      }
    }
    if(verbose){
      double absr=(alpha*mu).abs().maxCoeff();
      double maxe = max_err.maxCoeff();
      int asum=round(alpha.sum());
      double m_rel_li=rel_li(0);
      printf("%4d %+13.6e %1.9e %4d %1.9e %1.9e\n",(int)iter,lnZ(0),maxe,(int) asum,2.1,m_rel_li);
    }
    
    
    if(iter>itermax){
      printf("Maximum iteration number reached: %+0.2d \n",(int)iter);
      printf("The log variational lower bound of the last step increased by %+0.2e\n",(lnZ-lnZ0).maxCoeff());
      break;
    }
    iter=iter+1;
    // Rcout<<iter<<": worst error is "<<max_err.maxCoeff()<<" tolerance is "<<tolerance<<std::endl;
    // Rcout<<iter<<": size of max_err is: "<<max_err.size()<<std::endl;
    // Rcout<<iter<<": size of lnZ is: "<<lnZ.size()<<std::endl;
    // 
    
  }
  Rcpp::Rcout<<"mean lnZ is: "<<lnZ.mean()<<std::endl;
  return(Rcpp::DataFrame::create(_["logodds"]=logodds,_["sigb"]=sigma_beta,_["lnZ"]=lnZ.transpose()));
  
}


//[[Rcpp::export]]
Rcpp::List rss_varbvsr_alt_naive_grid(const c_sparseMatrix_internal SiRiS,
                                           const arrayxd_external sigma_beta,
                                           const arrayxd_external logodds,
                                           const arrayxd_external betahat,
                                           const arrayxd_external se,
                                           const arrayxd_external talpha0,
                                           const arrayxd_external tmu0,
                                           const arrayxd_external tSiRiSr0,
                                           double tolerance,
                                           int itermax,
                                           bool verbose,
                                           bool lnztol){
  
  
  //  bool lnztol=lnz_tol[0];
  using namespace Rcpp;
  using namespace Eigen;
  size_t sigb_size= sigma_beta.size();
  size_t logodds_size=logodds.size();
  size_t tot_size=sigb_size;
  if(tot_size!=logodds_size){
    Rcpp::stop("tot_size!=logodds_size");
  }
  size_t p=betahat.size();
  
  ArrayXXd alpha(p,tot_size);
  alpha.colwise()=talpha0;
  ArrayXXd mu(p,tot_size);
  mu.colwise()=tmu0;
  
  ArrayXXd SiRiSr(p,tot_size);
  SiRiSr.colwise()=tSiRiSr0;
  
  RowArray lnZ(tot_size);
  
  
  Eigen::ArrayXXd alpha0=alpha;
  Eigen::ArrayXXd mu0=mu;
//  Eigen::ArrayXXd SiRiSr0=SiRiSr;
  
  Eigen::ArrayXXd alpha1=alpha;
  Eigen::ArrayXXd mu1=mu;
//  Eigen::ArrayXXd SiRiSr1=SiRiSr;

  
  
  Eigen::ArrayXd sesquare =se*se;
  Eigen::ArrayXd  q= betahat/sesquare;
  Eigen::ArrayXXd s(p,tot_size);
  // Rcout<<"Initializing sigma_square and initial lnZ values"<<std::endl;
  for(size_t i=0; i<tot_size; i++){
    s.col(i)=(sesquare*(sigma_beta(i)*sigma_beta(i)))/(sesquare+(sigma_beta(i)*sigma_beta(i)));
    lnZ(i)=calculate_lnZ(q,alpha.col(i)*mu.col(i),SiRiSr.col(i),logodds(i),sesquare,alpha.col(i),mu.col(i),s.col(i),sigma_beta(i));
  }
  
  RowArray max_err(tot_size);
  max_err.setOnes();
  RowArray lnZ0=lnZ;
  RowArray rel_li(tot_size);
  rel_li.setZero();
  size_t iter=0;
  while(max_err.maxCoeff()>tolerance){
    //    Rcout<<iter<<": worst error is "<<max_err.maxCoeff()<<" tolerance is "<<tolerance<<std::endl;
    lnZ0=lnZ;
    alpha0=alpha;
    mu0=mu;
    
    bool reverse = iter%2!=0;
    
    rss_varbvsr_iter_alt_grid(SiRiS,sigma_beta,logodds,betahat,se,s,alpha,mu,SiRiSr,reverse);
    
    for(size_t c=0; c<tot_size; c++){
      lnZ(c)=calculate_lnZ(q,alpha.col(c)*mu.col(c),SiRiSr.col(c),logodds(c),sesquare,alpha.col(c),mu.col(c),s.col(c),sigma_beta(c));
    }
    //    rel_l0=rel_err(lnZ00,lnZ);
    rel_li=(lnZ-lnZ0).abs()/(lnZ.abs()+lnZ0.abs()+double_lim::epsilon());
    if(lnztol){
      max_err=rel_li;
    }else{
      for(size_t c=0; c<tot_size; c++){
        max_err(c)=find_maxerr(alpha.col(c),alpha0.col(c),alpha.col(c)*mu.col(c),alpha0.col(c)*mu0.col(c));
      }
    }
    iter++;
  }
  return Rcpp::List::create(Rcpp::Named("alpha")=alpha,
                            Rcpp::Named("mu")=mu,
                            Rcpp::Named("SiRiSr")=SiRiSr,
                            Rcpp::Named("max_err")=max_err,
                            Rcpp::Named("lnZ")=lnZ,
                            Rcpp::Named("iter")=iter);
}




//[[Rcpp::export]]
Rcpp::DataFrame grid_search_rss_varbvsr_alt_grid(const sparseMatrix_external SiRiS,
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
  
  return(rss_varbvsr_alt_grid(SiRiS,sigma_beta,logodds,betahat,se,talpha0,tmu0,tSiRiSr0,tolerance,itermax,verbose(0),lnz_tol(0)));
}

