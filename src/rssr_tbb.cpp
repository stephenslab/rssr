#include <RcppEigen.h>
#include "rssr.h"





Rcpp::NumericVector flatten_d(const fitdtype & obj,size_t retsize){

  Rcpp::NumericVector retvec(retsize);
  flatdtype dview= tbb::flatten2d(obj);
  size_t c=0;
  for(flatdtype::const_iterator i=dview.begin();i!=dview.end();i++){
    retvec[c]=*i;
    c++;
    if(c>retsize){
      Rcpp::Rcerr<<"Flattening failed, c:"<<c<<" retsize: "<<retsize<<std::endl;
      Rcpp::stop("flattening failed, too many items to flatten");
    }
  }
  if(c<=(retsize-1)){
    Rcpp::Rcerr<<"Flattening failed, c:"<<c<<" retsize: "<<retsize<<std::endl;
    Rcpp::stop("flattening failed, not enough items to flatten!");
  }
  return(retvec);
}
Rcpp::IntegerVector flatten_i(const fititype & obj,size_t retsize){

  Rcpp::IntegerVector retvec(retsize);
  flatitype dview= tbb::flatten2d(obj);
  size_t c=0;
  for(flatitype::const_iterator i=dview.begin();i!=dview.end();i++){
    retvec[c]=*i;
    c++;
    if(c>retsize){
      Rcpp::Rcerr<<"Flattening failed, c:"<<c<<" retsize: "<<retsize<<std::endl;
      Rcpp::stop("flattening failed, too many items to flatten");
    }
  }
  if(c<=(retsize-1)){
    Rcpp::Rcerr<<"Flattening failed, c:"<<c<<" retsize: "<<retsize<<std::endl;
    Rcpp::stop("flattening failed, not enough items to flatten!");
  }
  return(retvec);
}


//[[Rcpp::export]]
Rcpp::DataFrame grid_search_rss_varbvsr_alt(
    const  Matrix_external SiRiS,
    const arrayxd_external sigma_beta,
    const arrayxd_external logodds,
    const arrayxd_external  betahat,
    const arrayxd_external  se,
    const arrayxd_external talpha0,
    const arrayxd_external tmu0,
    const arrayxd_external tSiRiSr0,
    const double tolerance,
    const int itermax,
    Rcpp::LogicalVector lnz_tol,    const int n=1,const int grainsize=1){



  using namespace Rcpp;
  using namespace tbb;
  bool islnz_tol = lnz_tol(0);

  size_t p=betahat.size();
  size_t tot_size=logodds.size();
  

  Eigen::MatrixXd tSiRiS=SiRiS.eval();
  Eigen::ArrayXd tsigma_beta=sigma_beta.eval();
  Eigen::ArrayXd tsigma_beta_square=sigma_beta.eval();
  Eigen::ArrayXd tlogodds=logodds.eval();
  Eigen::ArrayXd tbetahat=betahat.eval();
  Eigen::ArrayXd tse=se.eval();

  

  Eigen::ArrayXd alpha0=talpha0.eval();

  

  
  
  Eigen::ArrayXd mu0=tmu0.eval();
  Eigen::ArrayXd SiRiSr0=tSiRiSr0.eval();



  double *all_alphasp,
    *all_alphas0p,
    *all_alphas1p,
    *all_musp,
    *all_mus0p,
    *all_mus1p,
    *all_SiRiSrp,
    *lnZp,
    *lnZ0p,
    *lnZ00p,
    *alpha_meanp,
    *mu_meanp,
    *final_max_errp,
    *mtpp,
    *pvevecp,
    *rvecp,
    *sesquarep,
    *qp,
    *sp,
    *ssratp;
  int* iternump,*final_btnump;
    

  all_alphasp = new double[p*tot_size];
  all_alphas0p = new double[p*tot_size];
  all_alphas1p = new double[p*tot_size];

  all_musp = new double[p*tot_size];
  all_mus0p = new double[p*tot_size];
  all_mus1p = new double[p*tot_size];

  all_SiRiSrp = new double[p*tot_size];
  lnZp = new double[tot_size];
  lnZ0p = new double[tot_size];
  lnZ00p = new double[tot_size];

  alpha_meanp = new double[tot_size];
  mu_meanp = new double[tot_size];
  final_max_errp = new double[tot_size];
  mtpp = new double[tot_size];
  pvevecp = new double[tot_size];
  rvecp = new double[tot_size];
  
  iternump = new int[tot_size];
  final_btnump = new int[tot_size];

  
  sesquarep = new double[p];
  qp = new double[p];
  sp = new double[tot_size*p];
  ssratp = new double[tot_size*p];
  
  
  

  mdarray sesquare(sesquarep,p);
  mdarray q(qp,p);
  m2darray s(sp,p,tot_size);
  m2darray ssrat(ssratp,p,tot_size);
  
  
  m2darray  all_alphas(all_alphasp,p,tot_size);
  m2darray  all_alphas0(all_alphas0p,p,tot_size);
  m2darray  all_alphas1(all_alphas1p,p,tot_size);
  
  m2darray  all_mus(all_musp,p,tot_size);
  m2darray  all_mus0(all_mus0p,p,tot_size);
  m2darray  all_mus1(all_mus1p,p,tot_size);
  m2darray  all_SiRiSr(all_SiRiSrp,p,tot_size);


  mdarray lnZ(lnZp,tot_size);
  mdarray lnZ0(lnZ0p,tot_size);
  mdarray lnZ00(lnZ00p,tot_size);
  
  mdarray alpha_mean(alpha_meanp,tot_size);
  mdarray mu_mean(mu_meanp,tot_size);
  mdarray final_max_err(final_max_errp,tot_size);
  miarray iternum(iternump,tot_size);
  miarray final_btnum(final_btnump,tot_size);

  mdarray pvevec(pvevecp,tot_size);
  mdarray rvec(rvecp,tot_size);
  mdarray mtp(mtpp,tot_size);
  
  
      
  
  for(size_t i=0; i<tot_size; i++){
    all_alphas.col(i)=talpha0;
    all_mus.col(i)=tmu0;
    all_SiRiSr.col(i)=tSiRiSr0;

  }
  

  Fit_wrap talpha(all_alphasp,all_alphas0p,all_alphas1p,p,tot_size);
  Fit_wrap tmu(all_musp,all_mus0p,all_mus1p,p,tot_size);
  Fit_wrap tSiRiSr(all_SiRiSrp,p,tot_size);

  
  Siris siris(tSiRiS);
  Data_wrap data_wrap(tbetahat.data(),tse.data(),sesquare.data(),q.data(),p);
  Param_wrap param(tlogodds.data(),tsigma_beta.data(),tsigma_beta_square.data(),s.data(),ssrat.data(),data_wrap,tot_size);
  Fit_res fit_res(lnZ.data(),lnZ0.data(),lnZ00.data(),alpha_mean.data(),
		  pvevec.data(),iternum.data(),mu_mean.data(),final_max_err.data(),rvec.data(),mtp.data(),final_btnum.data(),p,tot_size);
  
  std::vector<int> forward_range(p);


  
  for(size_t i=0;i<p;i++){
    forward_range[i]=i;
  }
  std::vector<int> reverse_range(p);
  std::reverse_copy(forward_range.begin(),forward_range.end(),reverse_range.begin());

  //  Rcpp::Rcout  <<"Initializing rssr_obj"<<std::endl;
    rssr rssr_obj(talpha,tmu,tSiRiSr,data_wrap,param,siris,fit_res,tolerance,n,itermax,forward_range,reverse_range);
    //  Rcpp::Rcout  <<"Starting parallel Algorithm!"<<std::endl;
  //rssr_obj(0,tot_size);
  //  Rcpp::Rcout<<"Starting Serial Algorithm!"<<std::endl;
  //  rssr_obj(blocked_range<size_t>(0,tot_size));
    parallel_for(blocked_range<size_t>(0,tot_size,grainsize),rssr_obj);
  //  Rcpp::Rcout<<"Finished!"<<std::endl;

  assert(logodds.size()==sigma_beta.size() &&" logodds==sigma_beta");
  assert(final_max_err.size()==sigma_beta.size() &&" logodds==sigma_beta");
  assert(iternum.size()==sigma_beta.size() &&" logodds==sigma_beta");

  NumericVector rlogodds(tot_size);
  NumericVector rsigma_beta(tot_size);
  NumericVector rfinal_max_err(tot_size);
  IntegerVector riternum(tot_size);
  NumericVector ralpha_mean(tot_size);

  NumericVector rmu_mean(tot_size);
  NumericVector rpvevec(tot_size);
  NumericVector rlnZ(tot_size);

  
  for(size_t i=0;i<tot_size;i++){
    rlogodds[i]=tlogodds(i);
    rsigma_beta[i]=tsigma_beta(i);
    rfinal_max_err[i]=final_max_err(i);
    riternum[i]=iternum(i);
    ralpha_mean[i]=alpha_mean(i);
    rmu_mean[i]=mu_mean(i);
    rpvevec[i]=pvevec(i);
    rlnZ[i]=lnZ(i);
  }

  
  // Rcpp::Rcout<<"logodds_sum:"<<logodds.sum()<<std::endl;
  // Rcpp::Rcout<<"sigma_beta_sum:"<<sigma_beta.sum()<<std::endl;
  // Rcpp::Rcout<<"iternum_sum:"<<iternum.sum()<<std::endl;
  // Rcpp::Rcout<<"alpha_mean_sum:"<<alpha_mean.sum()<<std::endl;
  // Rcpp::Rcout<<"mu_mean_sum:"<<mu_mean.sum()<<std::endl;
  // Rcpp::Rcout<<"pvevec_sum:"<<pvevec.sum()<<std::endl;
  // Rcpp::Rcout<<"lnZ_sum:"<<lnZ.sum()<<std::endl;
    
  //  Rcpp::Rcout<<"Returning"<<std::endl;
  assert(std::accumulate(rlogodds.begin(),rlogodds.end(),0.0)!=0 &&"rlogodds shouldn't sum to 0");

  assert(std::accumulate(rsigma_beta.begin(),rsigma_beta.end(),0.0)!=0 &&"rsigma_beta shouldn't sum to 0");
  assert(std::accumulate(rfinal_max_err.begin(),rfinal_max_err.end(),0.0)!=0 &&"rfinal_max_err shouldn't sum to 0");
  assert(std::accumulate(riternum.begin(),riternum.end(),0.0)!=0 &&"rfinal_max_err shouldn't sum to 0");
  assert(std::accumulate(ralpha_mean.begin(),ralpha_mean.end(),0.0)!=0 &&"ralpha_mean shouldn't sum to 0");
  assert(std::accumulate(rmu_mean.begin(),rmu_mean.end(),0.0)!=0 &&"rmu_mean shouldn't sum to 0");
  assert(std::accumulate(rmu_mean.begin(),rmu_mean.end(),0.0)!=0 &&"rmu_mean shouldn't sum to 0");
  assert(std::accumulate(rpvevec.begin(),rpvevec.end(),0.0)!=0 &&"rpvevec shouldn't sum to 0");
  assert(std::accumulate(rlnZ.begin(),rlnZ.end(),0.0)!=0 &&"rlnZ shouldn't sum to 0");
  
  
    
  
  
  
  
  delete [] all_alphasp;
  delete [] all_alphas0p;
  delete [] all_alphas1p;

  delete [] all_musp;
  delete [] all_mus0p;
  delete [] all_mus1p;

  delete [] all_SiRiSrp;
  delete [] lnZp;
  delete [] lnZ0p;
  delete [] lnZ00p;

  delete [] alpha_meanp;
  delete [] mu_meanp;
  delete [] final_max_errp;
  delete [] mtpp;
  delete [] pvevecp;
  delete [] rvecp ;
  
  delete [] iternump;
  delete [] final_btnump;

  
  delete []  sesquarep;
  delete []  qp;
  delete []  sp;
  delete []  ssratp;
  
  
    
  return(Rcpp::DataFrame::create(_["logodds"]=rlogodds,
				 _["sigb"]=rsigma_beta,
				 _["rel_err"]=rfinal_max_err,
				 _["iterations"]=riternum,
				 _["alpha_mean"]=ralpha_mean,
				 _["mu_mean"]=rmu_mean,
				 _["pve"]=rpvevec,
				 _["lnZ"]=rlnZ));
}






//[[Rcpp::export]]
Rcpp::DataFrame grid_search_rss_varbvsr_tls(
    const  Matrix_external SiRiS,
    const arrayxd_external sigma_beta,
    const arrayxd_external logodds,
    const arrayxd_external  betahat,
    const arrayxd_external  se,
    const arrayxd_external talpha0,
    const arrayxd_external tmu0,
    const arrayxd_external tSiRiSr0,
    const double tolerance,
    const int itermax,
    Rcpp::LogicalVector lnz_tol,
    const int n=1,const int grainsize=1){

  using namespace Rcpp;
  using namespace tbb;
  bool islnz_tol = lnz_tol(0);

  size_t p=betahat.size();
  size_t tot_size=logodds.size();
  

  Eigen::MatrixXd tSiRiS=SiRiS.eval();
  Eigen::ArrayXd tsigma_beta=sigma_beta.eval();
  Eigen::ArrayXd tsigma_beta_square=sigma_beta.eval();
  Eigen::ArrayXd tlogodds=logodds.eval();
  Eigen::ArrayXd tbetahat=betahat.eval();
  Eigen::ArrayXd tse=se.eval();
  Eigen::ArrayXd sesquare=se.eval();
  Eigen::ArrayXd q=se.eval();

  

  Eigen::ArrayXd alpha0=talpha0.eval();
  Eigen::ArrayXd mu0=tmu0.eval();
  Eigen::ArrayXd SiRiSr0=tSiRiSr0.eval();

  Fit_wrap talpha(alpha0.data(),alpha0.data(),alpha0.data(),p);
  Fit_wrap tmu(mu0.data(),mu0.data(),mu0.data(),p);
  Fit_wrap tSiRiSr(SiRiSr0.data(),p);

  
  Siris siris(tSiRiS.data(),p);
  Data_wrap data_wrap(tbetahat.data(),tse.data(),sesquare.data(),q.data(),p);
  Param_wrap param(tlogodds.data(),tsigma_beta.data(),tot_size,p);
  Fit_res fit_res(p,tot_size);
  
  std::vector<int> forward_range(p);
  for(size_t i=0;i<p;i++){
    forward_range[i]=i;
  }
  
  std::vector<int> reverse_range(p);
  std::reverse_copy(forward_range.begin(),forward_range.end(),reverse_range.begin());

  //  Rcpp::Rcout  <<"Initializing rssr_obj"<<std::endl;
  rssr rssr_obj(talpha,tmu,tSiRiSr,data_wrap,param,siris,fit_res,tolerance,n,itermax,forward_range,reverse_range);
  //Rcpp::Rcout  <<"Starting parallel Algorithm!"<<std::endl;
  //rssr_obj(0,tot_size);
  //  Rcpp::Rcout<<"Starting Serial Algorithm!"<<std::endl;
  //  rssr_obj(blocked_range<size_t>(0,tot_size));
  parallel_for(blocked_range<size_t>(0,tot_size,grainsize),rssr_obj);
  //    Rcpp::Rcout<<"Finished!"<<std::endl;
  
  assert(logodds.size()==sigma_beta.size() &&" logodds==sigma_beta");
  assert(final_max_err.size()==sigma_beta.size() &&" logodds==sigma_beta");
  assert(iternum.size()==sigma_beta.size() &&" logodds==sigma_beta");
  
  NumericVector rlogodds=flatten_d(rssr_obj.paraams.logodds_t,tot_size);
  NumericVector rsigma_beta=flatten_d(rssr_obj.paraams.sigma_beta_t,tot_size);
  NumericVector rfinal_max_err=flatten_d(rssr_obj.results.final_max_err_t,tot_size);
  IntegerVector riternum=flatten_i(rssr_obj.results.iternum_t,tot_size);
  NumericVector ralpha_mean=flatten_d(rssr_obj.results.alpha_mean_t,tot_size);
  NumericVector rmu_mean=flatten_d(rssr_obj.results.mu_mean_t,tot_size);
  NumericVector rpvevec=flatten_d(rssr_obj.results.pvevec_t,tot_size);
  NumericVector rlnZ=flatten_d(rssr_obj.results.lnZ_t,tot_size);


  
  assert(std::accumulate(rlogodds.begin(),rlogodds.end(),0.0)!=0 &&"rlogodds shouldn't sum to 0");

  assert(std::accumulate(rsigma_beta.begin(),rsigma_beta.end(),0.0)!=0 &&"rsigma_beta shouldn't sum to 0");
  assert(std::accumulate(rfinal_max_err.begin(),rfinal_max_err.end(),0.0)!=0 &&"rfinal_max_err shouldn't sum to 0");
  assert(std::accumulate(riternum.begin(),riternum.end(),0.0)!=0 &&"rfinal_max_err shouldn't sum to 0");
  assert(std::accumulate(ralpha_mean.begin(),ralpha_mean.end(),0.0)!=0 &&"ralpha_mean shouldn't sum to 0");
  assert(std::accumulate(rmu_mean.begin(),rmu_mean.end(),0.0)!=0 &&"rmu_mean shouldn't sum to 0");
  assert(std::accumulate(rmu_mean.begin(),rmu_mean.end(),0.0)!=0 &&"rmu_mean shouldn't sum to 0");
  assert(std::accumulate(rpvevec.begin(),rpvevec.end(),0.0)!=0 &&"rpvevec shouldn't sum to 0");
  assert(std::accumulate(rlnZ.begin(),rlnZ.end(),0.0)!=0 &&"rlnZ shouldn't sum to 0");

  
  
    
  return(Rcpp::DataFrame::create(_["logodds"]=rlogodds,
				 _["sigb"]=rsigma_beta,
				 _["rel_err"]=rfinal_max_err,
				 _["iterations"]=riternum,
				 _["alpha_mean"]=ralpha_mean,
				 _["mu_mean"]=rmu_mean,
				 _["pve"]=rpvevec,
				 _["lnZ"]=rlnZ));
}



				

//[[Rcpp::export]]
Rcpp::DataFrame grid_search_rss_varbvsr_norm_alt(
    const  Matrix_external SiRiS,
    const arrayxd_external sigma_beta,
    const arrayxd_external  betahat,
    const arrayxd_external  se,
    const arrayxd_external tmu0,
    const arrayxd_external tSiRiSr0,
    const double tolerance,
    const int itermax,
    Rcpp::LogicalVector lnz_tol,const int n=1,const int grainsize=1){



  using namespace Rcpp;
  using namespace tbb;
  bool islnz_tol = lnz_tol(0);

  size_t p=betahat.size();
  size_t tot_size=sigma_beta.size();
  

  Eigen::MatrixXd tSiRiS=SiRiS.eval();
  Eigen::ArrayXd tsigma_beta=sigma_beta.eval();
  Eigen::ArrayXd tsigma_beta_square=sigma_beta.eval();
  Eigen::ArrayXd tbetahat=betahat.eval();
  Eigen::ArrayXd tse=se.eval();
  Eigen::ArrayXd mu0=tmu0.eval();
  Eigen::ArrayXd SiRiSr0=tSiRiSr0.eval();



  double *all_musp,
    *all_mus0p,
    *all_mus1p,
    *all_SiRiSrp,
    *lnZp,
    *lnZ0p,
    *lnZ00p,
    *alpha_meanp,
    *mu_meanp,
    *final_max_errp,
    *mtpp,
    *pvevecp,
    *rvecp,
    *sesquarep,
    *qp,
    *sp,
    *ssratp;
  int* iternump,*final_btnump;

  all_musp = new double[p*tot_size];
  all_mus0p = new double[p*tot_size];
  all_mus1p = new double[p*tot_size];

  all_SiRiSrp = new double[p*tot_size];
  lnZp = new double[tot_size];
  lnZ0p = new double[tot_size];
  lnZ00p = new double[tot_size];

  alpha_meanp = new double[tot_size];
  mu_meanp = new double[tot_size];
  final_max_errp = new double[tot_size];
  mtpp = new double[tot_size];
  pvevecp = new double[tot_size];
  rvecp = new double[tot_size];
  
  iternump = new int[tot_size];
  final_btnump = new int[tot_size];

  
  sesquarep = new double[p];
  qp = new double[p];
  sp = new double[tot_size*p];
  ssratp = new double[tot_size*p];
  
  
  

  mdarray sesquare(sesquarep,p);
  mdarray q(qp,p);
  m2darray s(sp,p,tot_size);
  m2darray ssrat(ssratp,p,tot_size);
  
  
  m2darray  all_mus(all_musp,p,tot_size);
  m2darray  all_mus0(all_mus0p,p,tot_size);
  m2darray  all_mus1(all_mus1p,p,tot_size);
  m2darray  all_SiRiSr(all_SiRiSrp,p,tot_size);


  mdarray lnZ(lnZp,tot_size);
  mdarray lnZ0(lnZ0p,tot_size);
  mdarray lnZ00(lnZ00p,tot_size);
  
  mdarray alpha_mean(alpha_meanp,tot_size);
  mdarray mu_mean(mu_meanp,tot_size);
  mdarray final_max_err(final_max_errp,tot_size);
  miarray iternum(iternump,tot_size);
  miarray final_btnum(final_btnump,tot_size);

  mdarray pvevec(pvevecp,tot_size);
  mdarray rvec(rvecp,tot_size);
  mdarray mtp(mtpp,tot_size);
  
  
      
  
  for(size_t i=0; i<tot_size; i++){
    all_mus.col(i)=tmu0;
    all_SiRiSr.col(i)=tSiRiSr0;

  }
  
  Fit_wrap tmu(all_musp,all_mus0p,all_mus1p,p,tot_size);
  Fit_wrap tSiRiSr(all_SiRiSrp,p,tot_size);

  
  Siris siris(tSiRiS.data(),p);
  Data_wrap data_wrap(tbetahat.data(),tse.data(),sesquare.data(),q.data(),p);
  Param_wrap param(tsigma_beta.data(),tsigma_beta_square.data(),s.data(),ssrat.data(),data_wrap,tot_size);
  Fit_res fit_res(lnZ.data(),lnZ0.data(),lnZ00.data(),alpha_mean.data(),
		  pvevec.data(),iternum.data(),mu_mean.data(),final_max_err.data(),rvec.data(),mtp.data(),final_btnum.data(),p,tot_size);
  
  std::vector<int> forward_range(p);


  
  for(size_t i=0;i<p;i++){
    forward_range[i]=i;
  }
  std::vector<int> reverse_range(p);
  std::reverse_copy(forward_range.begin(),forward_range.end(),reverse_range.begin());

  Rcpp::Rcout  <<"Initializing rssr_obj"<<std::endl;
  rssr_norm rssr_obj(tmu,tSiRiSr,data_wrap,param,siris,fit_res,tolerance,n,itermax,forward_range,reverse_range);
  parallel_for(blocked_range<size_t>(0,tot_size,(size_t)tot_size/2),rssr_obj);



  assert(final_max_err.size()==sigma_beta.size() &&" logodds==sigma_beta");
  assert(iternum.size()==sigma_beta.size() &&" logodds==sigma_beta");

  NumericVector rsigma_beta(tot_size);
  NumericVector rfinal_max_err(tot_size);
  IntegerVector riternum(tot_size);
  NumericVector ralpha_mean(tot_size);

  NumericVector rmu_mean(tot_size);
  NumericVector rpvevec(tot_size);
  NumericVector rlnZ(tot_size);

  
  for(size_t i=0;i<tot_size;i++){
    rsigma_beta[i]=tsigma_beta(i);
    rfinal_max_err[i]=final_max_err(i);
    riternum[i]=iternum(i);
    ralpha_mean[i]=alpha_mean(i);
    rmu_mean[i]=mu_mean(i);
    rpvevec[i]=pvevec(i);
    rlnZ[i]=lnZ(i);
  }



  assert(std::accumulate(rsigma_beta.begin(),rsigma_beta.end(),0.0)!=0 &&"rsigma_beta shouldn't sum to 0");
  assert(std::accumulate(rfinal_max_err.begin(),rfinal_max_err.end(),0.0)!=0 &&"rfinal_max_err shouldn't sum to 0");
  assert(std::accumulate(riternum.begin(),riternum.end(),0.0)!=0 &&"rfinal_max_err shouldn't sum to 0");
  assert(std::accumulate(ralpha_mean.begin(),ralpha_mean.end(),0.0)!=0 &&"ralpha_mean shouldn't sum to 0");
  assert(std::accumulate(rmu_mean.begin(),rmu_mean.end(),0.0)!=0 &&"rmu_mean shouldn't sum to 0");
  assert(std::accumulate(rmu_mean.begin(),rmu_mean.end(),0.0)!=0 &&"rmu_mean shouldn't sum to 0");
  assert(std::accumulate(rpvevec.begin(),rpvevec.end(),0.0)!=0 &&"rpvevec shouldn't sum to 0");
  assert(std::accumulate(rlnZ.begin(),rlnZ.end(),0.0)!=0 &&"rlnZ shouldn't sum to 0");
  
  

  delete [] all_musp;
  delete [] all_mus0p;
  delete [] all_mus1p;

  delete [] all_SiRiSrp;
  delete [] lnZp;
  delete [] lnZ0p;
  delete [] lnZ00p;

  delete [] alpha_meanp;
  delete [] mu_meanp;
  delete [] final_max_errp;
  delete [] mtpp;
  delete [] pvevecp;
  delete [] rvecp ;
  
  delete [] iternump;
  delete [] final_btnump;

  
  delete []  sesquarep;
  delete []  qp;
  delete []  sp;
  delete []  ssratp;
  
  
    
  return(Rcpp::DataFrame::create(_["sigb"]=rsigma_beta,
				 _["rel_err"]=rfinal_max_err,
				 _["iterations"]=riternum,
				 _["alpha_mean"]=ralpha_mean,
				 _["mu_mean"]=rmu_mean,
				 _["pve"]=rpvevec,
				 _["lnZ"]=rlnZ));
}
