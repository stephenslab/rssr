#include <RcppEigen.h>
#include "rssr.h"


//# define rassert(EX) (void)((EX) || (pkassert (#EX, __FILE__, __LINE__),0))






Rcpp::NumericMatrix initialize_smat(const Rcpp::NumericVector &sesquare, const Rcpp::NumericVector &sigma_beta_square){
  c_mdarray tsesq(&sesquare[0],(size_t)sesquare.size());

  Rcpp::NumericMatrix retmat(tsesq.size(),sigma_beta_square.size());
  mmat tretmat(&retmat[0],retmat.nrow(),retmat.ncol());

  for(size_t i=0; i<retmat.ncol();i++) {
    tretmat.col(i)=(tsesq*(sigma_beta_square[i]))/(tsesq+(sigma_beta_square[i]));
  }
  rassert(isnan(tretmat.array()).sum()==0 && "s has no NaN  values");
  rassert(tretmat.sum()!=0 && "smat shouldn't be 0");
  return(retmat);
}


Rcpp::NumericMatrix initialize_ssratmat(const Rcpp::NumericMatrix &s, const Rcpp::NumericVector &sigma_beta_square){
  c_mmat ts(&s[0],s.nrow(),s.ncol());
  c_mdarray tsigsq(&sigma_beta_square[0],sigma_beta_square.size());
  Rcpp::NumericMatrix retmat(s.nrow(),s.ncol());
  mmat tretmat(&retmat[0],retmat.nrow(),retmat.ncol());
  rassert(isnan(ts.array()).sum()==0 && "s has no NaN  values");
  rassert(tsigsq.minCoeff()>0 && "sigma_beta_square has all non-negative values");
  rassert(isnan(tsigsq.array()).sum()==0 && "sigma_beta_square has no NaN  values");
  
  for(size_t i=0; i<s.ncol();i++){
    tretmat.col(i)=(ts.col(i).array()/tsigsq(i)).log();
  }
  
  rassert(isnan(tretmat.array()).sum()==0 && "ssrat has no NaN  values");

  return(retmat);
}



//[[Rcpp::export(name="grid_search_rss_varbvsr_tls")]]
Rcpp::DataFrame grid_search_rss_varbvsr_tls_exp(
    const  Rcpp::NumericMatrix &SiRiS,
    const Rcpp::NumericVector &sigma_beta,
    const Rcpp::NumericVector &logodds,
    const Rcpp::NumericVector &betahat,
    const Rcpp::NumericVector &se,
    const Rcpp::NumericVector &alpha0,
    const Rcpp::NumericVector &mu0,
    const Rcpp::NumericVector &SiRiSr0,
    const double tolerance,
    const int itermax,
    Rcpp::LogicalVector lnz_tol,
    const int n=1,const int grainsize=1){

  return(grid_search_rss_varbvsr_tls(SiRiS,
				     sigma_beta,
				     logodds,
				     betahat,
				     se,
				     alpha0,
				     mu0,
				     SiRiSr0,
				     tolerance,
				     itermax,
				     lnz_tol,
				     n, grainsize));
}


Rcpp::DataFrame grid_search_rss_varbvsr_tls(
    const  Rcpp::NumericMatrix &SiRiS,
    const Rcpp::NumericVector &sigma_beta,
    const Rcpp::NumericVector &logodds,
    const Rcpp::NumericVector &betahat,
    const Rcpp::NumericVector &se,
    const Rcpp::NumericVector &alpha0,
    const Rcpp::NumericVector &mu0,
    const Rcpp::NumericVector &SiRiSr0,
    const double tolerance,
    const int itermax,
    Rcpp::LogicalVector lnz_tol,
    const int n=1,const int grainsize=1){

  using namespace Rcpp;
  using namespace tbb;
  using namespace RcppParallel;
  bool islnz_tol = lnz_tol(0);

  size_t p=betahat.size();
  size_t tot_size=logodds.size();


  
  NumericVector sigma_beta_square(tot_size);
  NumericVector sesquare(p);
  NumericVector q(p);
  
  std::transform(sigma_beta.begin(),sigma_beta.end(),sigma_beta_square.begin(),[](double c) {return c*c;});
  std::transform(se.begin(),se.end(),sesquare.begin(),[](double c) {return c*c;});
  std::transform(betahat.begin(),betahat.end(),sesquare.begin(),q.begin(),std::divides<double>());


  NumericMatrix s=initialize_smat(sesquare,sigma_beta_square);
  NumericMatrix ssrat=initialize_ssratmat(s,sigma_beta_square);
  

  RMatrix<double> tSiRiS(SiRiS);
  RVector<double> tsigma_beta(sigma_beta);
  RVector<double> tsigma_beta_square(sigma_beta_square);
  RVector<double> tlogodds(logodds);
  RVector<double> tbetahat(betahat);
  RVector<double> tse(se);
  RVector<double> tsesquare(sesquare);
  RVector<double> tq(q);

  RMatrix<double> ts(s);
  RMatrix<double> tssrat(ssrat);


  RVector<double> talpha0(alpha0);
  RVector<double> tmu0(mu0);
  RVector<double> tSiRiSr0(tSiRiSr0);

  NumericMatrix alphas(p,tot_size);
  RMatrix<double>talphas(alphas);
  
  NumericMatrix alpha0s(p,tot_size);
  RMatrix<double>talpha0s(alpha0s);
  
  NumericMatrix alpha1s(p,tot_size);
  RMatrix<double>talpha1s(alpha1s);
  
  NumericMatrix mus(p,tot_size);
  RMatrix<double>tmus(mus);
  
  NumericMatrix mu0s(p,tot_size);
  RMatrix<double>tmu0s(mu0s);
  
  NumericMatrix mu1s(p,tot_size);
  RMatrix<double>tmu1s(mu1s);
  

  NumericMatrix SiRiSrs(p,tot_size);
  RMatrix<double>tSiRiSrs(SiRiSrs);
  

  mdarray ttalpha0(talpha0.begin(),p);
  rassert(ttalpha0.minCoeff()>0 && "alpha0 has all non-negative values");

  Fit_wrap talpha(alpha0.begin(),alphas.begin(),alpha0s.begin(),alpha1s.begin(),p,tot_size);
  Fit_wrap tmu(mu0.begin(),mus.begin(),mu0s.begin(),mu1s.begin(),p,tot_size);
  Fit_wrap tSiRiSr(SiRiSr0.begin(),tSiRiSrs.begin(),p,tot_size);
  
  m2darray ttalpha(talpha.fit_p,p,tot_size);
  //  std::cout<<"ttalph.amincoeff: "<<ttalpha.minCoeff()<<std::endl;
  //  std::cout<<"ttalpha0.mincoeff: "<<ttalpha0.minCoeff()<<std::endl;
  rassert(ttalpha.minCoeff()>0 && "alpha has all non-negative values");


  // NumericVector rlogodds(tot_size);
  // NumericVector rsigma_beta(tot_size);
  NumericVector rfinal_max_err(tot_size);
  IntegerVector rfinal_btnum(tot_size);
  IntegerVector riternum(tot_size);
  NumericVector ralpha_mean(tot_size);
  NumericVector rmu_mean(tot_size);
  NumericVector rpvevec(tot_size);
  NumericVector rlnZ(tot_size);

  // std::fill(lnZ.begin(),lnZ.end(),0);
  // std::fill(lnZ0.begin(),lnZ.end(),0);

  NumericVector lnZ0(tot_size);
  
  NumericVector lnZ00(tot_size);

  NumericVector rvec(tot_size);
  NumericVector mtp(tot_size);
  
  RVector<double> tlnZ0(lnZ0);
  RVector<double> tlnZ00(lnZ00);
  RVector<double> trvec(rvec);
  RVector<double> tmtp(mtp);
  

  // RVector<double> trlogodds(rlogodds);
  // RVector<double> trsigma_beta(rsigma_beta);
  RVector<double> trfinal_max_err(rfinal_max_err);
  RVector<int> trfinal_btnum(rfinal_btnum);
  IntegerVector triternum(riternum);
  RVector<double> tralpha_mean(ralpha_mean);
  RVector<double> trmu_mean(rmu_mean);
  RVector<double> trpvevec(rpvevec);
  RVector<double> trlnZ(rlnZ);
  

  
  Data_wrap data_wrap(tbetahat.begin(),tse.begin(),sesquare.begin(),q.begin(),tSiRiS.begin(),p);
  Param_wrap param(tlogodds.begin(),tsigma_beta.begin(),tsigma_beta_square.begin(),ts.begin(),tssrat.begin(),tot_size,p);
  Fit_res fit_res(trlnZ.begin(),
		  tlnZ0.begin(),
		  tlnZ00.begin(),
		  tralpha_mean.begin(),
		  trmu_mean.begin(),
		  trpvevec.begin(),
		  triternum.begin(),
		  trfinal_max_err.begin(),
		  trfinal_btnum.begin(),
		  trvec.begin(),tmtp.begin(),
		  p,
		  tot_size);
  
  std::vector<int> forward_range(p);
  for(size_t i=0;i<p;i++){
    forward_range[i]=i;
  }
  
  std::vector<int> reverse_range(p);
  std::reverse_copy(forward_range.begin(),forward_range.end(),reverse_range.begin());

  //  Rcpp::Rcout  <<"Initializing rssr_obj"<<std::endl;
  rssr rssr_obj(&talpha,&tmu,&tSiRiSr,&data_wrap,&param,&fit_res,tolerance,n,itermax,forward_range,reverse_range);
  //Rcpp::Rcout  <<"Starting parallel Algorithm!"<<std::endl;
  //rssr_obj(0,tot_size);
  //Rcpp::Rcout<<"Starting Serial Algorithm!"<<std::endl;
  //  rssr_obj(blocked_range<size_t>(0,tot_size));
  parallel_for(blocked_range<size_t>(0,tot_size,grainsize),rssr_obj);
  //    Rcpp::Rcout<<"Finished!"<<std::endl;
  //  std::cout<<"Checking then Returning!"<<std::endl;
  assert(logodds.size()==sigma_beta.size() &&" logodds==sigma_beta");
  assert(final_max_err.size()==sigma_beta.size() &&" logodds==sigma_beta");
  assert(iternum.size()==sigma_beta.size() &&" logodds==sigma_beta");
  


  


  
  assert(std::accumulate(rlogodds.begin(),rlogodds.end(),0.0)!=0 &&"rlogodds shouldn't sum to 0");

  assert(std::accumulate(rsigma_beta.begin(),rsigma_beta.end(),0.0)!=0 &&"rsigma_beta shouldn't sum to 0");
  assert(std::accumulate(rfinal_max_err.begin(),rfinal_max_err.end(),0.0)!=0 &&"rfinal_max_err shouldn't sum to 0");
  assert(std::accumulate(riternum.begin(),riternum.end(),0.0)!=0 &&"rfinal_max_err shouldn't sum to 0");
  assert(std::accumulate(ralpha_mean.begin(),ralpha_mean.end(),0.0)!=0 &&"ralpha_mean shouldn't sum to 0");
  assert(std::accumulate(rmu_mean.begin(),rmu_mean.end(),0.0)!=0 &&"rmu_mean shouldn't sum to 0");
  assert(std::accumulate(rmu_mean.begin(),rmu_mean.end(),0.0)!=0 &&"rmu_mean shouldn't sum to 0");
  assert(std::accumulate(rpvevec.begin(),rpvevec.end(),0.0)!=0 &&"rpvevec shouldn't sum to 0");
  assert(std::accumulate(rlnZ.begin(),rlnZ.end(),0.0)!=0 &&"rlnZ shouldn't sum to 0");

  //    std::cout<<"Returning!"<<std::endl;
  
    
  return(Rcpp::DataFrame::create(_["logodds"]=logodds,
				 _["sigb"]=sigma_beta,
				 _["rel_err"]=rfinal_max_err,
				 _["iterations"]=riternum,
				 _["alpha_mean"]=ralpha_mean,
				 _["mu_mean"]=rmu_mean,
				 _["pve"]=rpvevec,
				 _["lnZ"]=rlnZ));
}



//[[Rcpp::export(name="grid_search_rss_varbvsr_norm_tls")]]
Rcpp::DataFrame grid_search_rss_varbvsr_norm_tls_exp(
    const  Rcpp::NumericMatrix &SiRiS,
    const Rcpp::NumericVector &sigma_beta,
    const Rcpp::NumericVector &betahat,
    const Rcpp::NumericVector &se,
    const Rcpp::NumericVector &mu0,
    const Rcpp::NumericVector &SiRiSr0,
    const double tolerance,
    const int itermax,
    Rcpp::LogicalVector lnz_tol,
    const int n=1,const int grainsize=1){
  
  return(grid_search_rss_varbvsr_norm_tls(SiRiS,
                                          sigma_beta,
                                          betahat,
                                          se,
                                          mu0,
                                          SiRiSr0,
                                          tolerance,
                                          itermax,
                                          lnz_tol,
                                          n, grainsize));
}


Rcpp::DataFrame grid_search_rss_varbvsr_norm_tls(const  Rcpp::NumericMatrix &SiRiS,
                                                 const Rcpp::NumericVector &sigma_beta,
                                                 const Rcpp::NumericVector &betahat,
                                                 const Rcpp::NumericVector &se,
                                                 const Rcpp::NumericVector &mu0,
                                                 const Rcpp::NumericVector &SiRiSr0,
                                                 const double tolerance,
                                                 const int itermax,
                                                 Rcpp::LogicalVector lnz_tol,
                                                 const int n=1,const int grainsize=1){
  
  
  using namespace Rcpp;
  using namespace tbb;
  using namespace RcppParallel;
  bool islnz_tol = lnz_tol(0);
  
  size_t p=betahat.size();
  size_t tot_size=sigma_beta.size();
  
  
  
  NumericVector sigma_beta_square(tot_size);
  NumericVector sesquare(p);
  NumericVector q(p);
  
  std::transform(sigma_beta.begin(),sigma_beta.end(),sigma_beta_square.begin(),[](double c) {return c*c;});
  std::transform(se.begin(),se.end(),sesquare.begin(),[](double c) {return c*c;});
  std::transform(betahat.begin(),betahat.end(),sesquare.begin(),q.begin(),std::divides<double>());
  
  
  NumericMatrix s=initialize_smat(sesquare,sigma_beta_square);
  NumericMatrix ssrat=initialize_ssratmat(s,sigma_beta_square);
  
  
  RMatrix<double> tSiRiS(SiRiS);
  RVector<double> tsigma_beta(sigma_beta);
  RVector<double> tsigma_beta_square(sigma_beta_square);
  RVector<double> tbetahat(betahat);
  RVector<double> tse(se);
  RVector<double> tsesquare(sesquare);
  RVector<double> tq(q);
  
  RMatrix<double> ts(s);
  RMatrix<double> tssrat(ssrat);
  
  
  
  RVector<double> tmu0(mu0);
  RVector<double> tSiRiSr0(SiRiSr0);
  
  NumericMatrix mus(p,tot_size);
  RMatrix<double>tmus(mus);
  
  NumericMatrix mu0s(p,tot_size);
  RMatrix<double>tmu0s(mu0s);
  
  NumericMatrix mu1s(p,tot_size);
  RMatrix<double>tmu1s(mu1s);
  
  
  NumericMatrix SiRiSrs(p,tot_size);
  RMatrix<double>tSiRiSrs(SiRiSrs);
  
  Fit_wrap tmu(mu0.begin(),mus.begin(),mu0s.begin(),mu1s.begin(),p,tot_size);
  Fit_wrap tSiRiSr(SiRiSr0.begin(),tSiRiSrs.begin(),p,tot_size);
  
  // NumericVector rlogodds(tot_size);
  // NumericVector rsigma_beta(tot_size);
  NumericVector rfinal_max_err(tot_size);
  IntegerVector rfinal_btnum(tot_size);
  IntegerVector riternum(tot_size);
  NumericVector ralpha_mean(tot_size);
  NumericVector rmu_mean(tot_size);
  NumericVector rpvevec(tot_size);
  NumericVector rlnZ(tot_size);
  
  
  NumericVector lnZ0(tot_size);
  NumericVector lnZ00(tot_size);
  
  NumericVector rvec(tot_size);
  NumericVector mtp(tot_size);
  
  RVector<double> tlnZ0(lnZ0);
  RVector<double> tlnZ00(lnZ00);
  RVector<double> trvec(rvec);
  RVector<double> tmtp(mtp);
  
  
  RVector<double> trfinal_max_err(rfinal_max_err);
  RVector<int> trfinal_btnum(rfinal_btnum);
  IntegerVector triternum(riternum);
  RVector<double> tralpha_mean(ralpha_mean);
  RVector<double> trmu_mean(rmu_mean);
  RVector<double> trpvevec(rpvevec);
  RVector<double> trlnZ(rlnZ);
  
  
  
  Data_wrap data_wrap(tbetahat.begin(),tse.begin(),sesquare.begin(),q.begin(),tSiRiS.begin(),p);
  Param_wrap param(tsigma_beta.begin(),tsigma_beta_square.begin(),ts.begin(),tssrat.begin(),tot_size,p);
  Fit_res fit_res(trlnZ.begin(),
                  tlnZ0.begin(),
                  tlnZ00.begin(),
                  tralpha_mean.begin(),
                  trmu_mean.begin(),
                  trpvevec.begin(),
                  triternum.begin(),
                  trfinal_max_err.begin(),
                  trfinal_btnum.begin(),
                  trvec.begin(),tmtp.begin(),
                  p,
                  tot_size);
  
  std::vector<int> forward_range(p);
  for(size_t i=0;i<p;i++){
    forward_range[i]=i;
  }
  
  std::vector<int> reverse_range(p);
  std::reverse_copy(forward_range.begin(),forward_range.end(),reverse_range.begin());
  
  //  Rcpp::Rcout  <<"Initializing rssr_obj"<<std::endl;
  rssr_norm rssr_obj(&tmu,&tSiRiSr,&data_wrap,&param,&fit_res,tolerance,n,itermax,forward_range,reverse_range);
  //Rcpp::Rcout  <<"Starting parallel Algorithm!"<<std::endl;
  //rssr_obj(0,tot_size);
  //Rcpp::Rcout<<"Starting Serial Algorithm!"<<std::endl;
  //  rssr_obj(blocked_range<size_t>(0,tot_size));
  parallel_for(blocked_range<size_t>(0,tot_size,grainsize),rssr_obj);
  //    Rcpp::Rcout<<"Finished!"<<std::endl;
  //  std::cout<<"Checking then Returning!"<<std::endl;
  assert(final_max_err.size()==sigma_beta.size() &&" final_max_err.size()==sigma_beta.size()");
  assert(iternum.size()==sigma_beta.size() &&" iternum.size()==sigma_beta.size()");
  
  
  
  
  
  
  
  
  assert(std::accumulate(rsigma_beta.begin(),rsigma_beta.end(),0.0)!=0 &&"rsigma_beta shouldn't sum to 0");
  assert(std::accumulate(rfinal_max_err.begin(),rfinal_max_err.end(),0.0)!=0 &&"rfinal_max_err shouldn't sum to 0");
  assert(std::accumulate(riternum.begin(),riternum.end(),0.0)!=0 &&"rfinal_max_err shouldn't sum to 0");
  assert(std::accumulate(ralpha_mean.begin(),ralpha_mean.end(),0.0)!=0 &&"ralpha_mean shouldn't sum to 0");
  assert(std::accumulate(rmu_mean.begin(),rmu_mean.end(),0.0)!=0 &&"rmu_mean shouldn't sum to 0");
  assert(std::accumulate(rmu_mean.begin(),rmu_mean.end(),0.0)!=0 &&"rmu_mean shouldn't sum to 0");
  assert(std::accumulate(rpvevec.begin(),rpvevec.end(),0.0)!=0 &&"rpvevec shouldn't sum to 0");
  assert(std::accumulate(rlnZ.begin(),rlnZ.end(),0.0)!=0 &&"rlnZ shouldn't sum to 0");
  
  //    std::cout<<"Returning!"<<std::endl;
  
  
  return(Rcpp::DataFrame::create(_["sigb"]=sigma_beta,
                                 _["rel_err"]=rfinal_max_err,
                                 _["iterations"]=riternum,
                                 _["alpha_mean"]=ralpha_mean,
                                 _["mu_mean"]=rmu_mean,
                                 _["pve"]=rpvevec,
                                 _["lnZ"]=rlnZ));
  
  
}



Rcpp::DataFrame grid_search_rss_varbvsr_norm_tls_sparse(const  Eigen::Map< SparseMatrix<double> >SiRiS,
                                                 const Rcpp::NumericVector &sigma_beta,
                                                 const Rcpp::NumericVector &betahat,
                                                 const Rcpp::NumericVector &se,
                                                 const Rcpp::NumericVector &mu0,
                                                 const Rcpp::NumericVector &SiRiSr0,
                                                 const double tolerance,
                                                 const int itermax,
                                                 Rcpp::LogicalVector lnz_tol,
                                                 const int n=1,const int grainsize=1){
  
  
  using namespace Rcpp;
  using namespace tbb;
  using namespace RcppParallel;
  bool islnz_tol = lnz_tol(0);
  
  size_t p=betahat.size();
  size_t tot_size=sigma_beta.size();
  
  
  
  NumericVector sigma_beta_square(tot_size);
  NumericVector sesquare(p);
  NumericVector q(p);
  
  std::transform(sigma_beta.begin(),sigma_beta.end(),sigma_beta_square.begin(),[](double c) {return c*c;});
  std::transform(se.begin(),se.end(),sesquare.begin(),[](double c) {return c*c;});
  std::transform(betahat.begin(),betahat.end(),sesquare.begin(),q.begin(),std::divides<double>());
  
  
  NumericMatrix s=initialize_smat(sesquare,sigma_beta_square);
  NumericMatrix ssrat=initialize_ssratmat(s,sigma_beta_square);
  
  
//  RMatrix<double> tSiRiS(SiRiS);
  RVector<double> tsigma_beta(sigma_beta);
  RVector<double> tsigma_beta_square(sigma_beta_square);
  RVector<double> tbetahat(betahat);
  RVector<double> tse(se);
  RVector<double> tsesquare(sesquare);
  RVector<double> tq(q);
  
  RMatrix<double> ts(s);
  RMatrix<double> tssrat(ssrat);
  
  
  
  RVector<double> tmu0(mu0);
  RVector<double> tSiRiSr0(SiRiSr0);
  
  NumericMatrix mus(p,tot_size);
  RMatrix<double>tmus(mus);
  
  NumericMatrix mu0s(p,tot_size);
  RMatrix<double>tmu0s(mu0s);
  
  NumericMatrix mu1s(p,tot_size);
  RMatrix<double>tmu1s(mu1s);
  
  
  NumericMatrix SiRiSrs(p,tot_size);
  RMatrix<double>tSiRiSrs(SiRiSrs);
  
  Fit_wrap tmu(mu0.begin(),mus.begin(),mu0s.begin(),mu1s.begin(),p,tot_size);
  Fit_wrap tSiRiSr(SiRiSr0.begin(),tSiRiSrs.begin(),p,tot_size);
  
  // NumericVector rlogodds(tot_size);
  // NumericVector rsigma_beta(tot_size);
  NumericVector rfinal_max_err(tot_size);
  IntegerVector rfinal_btnum(tot_size);
  IntegerVector riternum(tot_size);
  NumericVector ralpha_mean(tot_size);
  NumericVector rmu_mean(tot_size);
  NumericVector rpvevec(tot_size);
  NumericVector rlnZ(tot_size);
  
  
  NumericVector lnZ0(tot_size);
  NumericVector lnZ00(tot_size);
  
  NumericVector rvec(tot_size);
  NumericVector mtp(tot_size);
  
  RVector<double> tlnZ0(lnZ0);
  RVector<double> tlnZ00(lnZ00);
  RVector<double> trvec(rvec);
  RVector<double> tmtp(mtp);
  
  
  RVector<double> trfinal_max_err(rfinal_max_err);
  RVector<int> trfinal_btnum(rfinal_btnum);
  IntegerVector triternum(riternum);
  RVector<double> tralpha_mean(ralpha_mean);
  RVector<double> trmu_mean(rmu_mean);
  RVector<double> trpvevec(rpvevec);
  RVector<double> trlnZ(rlnZ);
  
  
  
  Data_wrap data_wrap(tbetahat.begin(),tse.begin(),sesquare.begin(),q.begin(),SiRiS,p);
  Param_wrap param(tsigma_beta.begin(),tsigma_beta_square.begin(),ts.begin(),tssrat.begin(),tot_size,p);
  Fit_res fit_res(trlnZ.begin(),
                  tlnZ0.begin(),
                  tlnZ00.begin(),
                  tralpha_mean.begin(),
                  trmu_mean.begin(),
                  trpvevec.begin(),
                  triternum.begin(),
                  trfinal_max_err.begin(),
                  trfinal_btnum.begin(),
                  trvec.begin(),tmtp.begin(),
                  p,
                  tot_size);
  
  std::vector<int> forward_range(p);
  for(size_t i=0;i<p;i++){
    forward_range[i]=i;
  }
  
  std::vector<int> reverse_range(p);
  std::reverse_copy(forward_range.begin(),forward_range.end(),reverse_range.begin());
  
  //  Rcpp::Rcout  <<"Initializing rssr_obj"<<std::endl;
  rssr_norm rssr_obj(&tmu,&tSiRiSr,&data_wrap,&param,&fit_res,tolerance,n,itermax,forward_range,reverse_range);
  //Rcpp::Rcout  <<"Starting parallel Algorithm!"<<std::endl;
  //rssr_obj(0,tot_size);
  //Rcpp::Rcout<<"Starting Serial Algorithm!"<<std::endl;
  //  rssr_obj(blocked_range<size_t>(0,tot_size));
  parallel_for(blocked_range<size_t>(0,tot_size,grainsize),rssr_obj);
  //    Rcpp::Rcout<<"Finished!"<<std::endl;
  //  std::cout<<"Checking then Returning!"<<std::endl;
  assert(final_max_err.size()==sigma_beta.size() &&" final_max_err.size()==sigma_beta.size()");
  assert(iternum.size()==sigma_beta.size() &&" iternum.size()==sigma_beta.size()");
  
  
  
  
  
  
  
  
  assert(std::accumulate(rsigma_beta.begin(),rsigma_beta.end(),0.0)!=0 &&"rsigma_beta shouldn't sum to 0");
  assert(std::accumulate(rfinal_max_err.begin(),rfinal_max_err.end(),0.0)!=0 &&"rfinal_max_err shouldn't sum to 0");
  assert(std::accumulate(riternum.begin(),riternum.end(),0.0)!=0 &&"rfinal_max_err shouldn't sum to 0");
  assert(std::accumulate(ralpha_mean.begin(),ralpha_mean.end(),0.0)!=0 &&"ralpha_mean shouldn't sum to 0");
  assert(std::accumulate(rmu_mean.begin(),rmu_mean.end(),0.0)!=0 &&"rmu_mean shouldn't sum to 0");
  assert(std::accumulate(rmu_mean.begin(),rmu_mean.end(),0.0)!=0 &&"rmu_mean shouldn't sum to 0");
  assert(std::accumulate(rpvevec.begin(),rpvevec.end(),0.0)!=0 &&"rpvevec shouldn't sum to 0");
  assert(std::accumulate(rlnZ.begin(),rlnZ.end(),0.0)!=0 &&"rlnZ shouldn't sum to 0");
  
  //    std::cout<<"Returning!"<<std::endl;
  
  
  return(Rcpp::DataFrame::create(_["sigb"]=sigma_beta,
                                 _["rel_err"]=rfinal_max_err,
                                 _["iterations"]=riternum,
                                 _["alpha_mean"]=ralpha_mean,
                                 _["mu_mean"]=rmu_mean,
                                 _["pve"]=rpvevec,
                                 _["lnZ"]=rlnZ));
  
  
}



//[[Rcpp::export(name="grid_search_rss_varbvsr_norm_sparse")]]
Rcpp::DataFrame grid_search_rss_varbvsr_norm_sparse_exp(
    const  Eigen::Map<Eigen::SparseMatrix<double> > SiRiS,
    const Rcpp::NumericVector &sigma_beta,
    const Rcpp::NumericVector &betahat,
    const Rcpp::NumericVector &se,
    const Rcpp::NumericVector &mu0,
    const Rcpp::NumericVector &SiRiSr0,
    const double tolerance,
    const int itermax,
    Rcpp::LogicalVector lnz_tol,
    const int n=1,const int grainsize=1){
  
  return(grid_search_rss_varbvsr_norm_tls_sparse(SiRiS,
                                          sigma_beta,
                                          betahat,
                                          se,
                                          mu0,
                                          SiRiSr0,
                                          tolerance,
                                          itermax,
                                          lnz_tol,
                                          n, grainsize));
}




//[[Rcpp::export(name="grid_search_rssr_varbvsr_sparse")]]
Rcpp::DataFrame grid_search_rss_varbvsr_sparse_exp(
    const  Eigen::Map<Eigen::SparseMatrix<double> > SiRiS,
    const Rcpp::NumericVector &sigma_beta,
    const Rcpp::NumericVector &logodds,
    const Rcpp::NumericVector &betahat,
    const Rcpp::NumericVector &se,
    const Rcpp::NumericVector &alpha0,
    const Rcpp::NumericVector &mu0,
    const Rcpp::NumericVector &SiRiSr0,
    const double tolerance,
    const int itermax,
    Rcpp::LogicalVector lnz_tol,
    const int n=1,const int grainsize=1){
  
  return(grid_search_rss_varbvsr_sparse(SiRiS,
                                     sigma_beta,
                                     logodds,
                                     betahat,
                                     se,
                                     alpha0,
                                     mu0,
                                     SiRiSr0,
                                     tolerance,
                                     itermax,
                                     lnz_tol,
                                     n, grainsize));
}


Rcpp::DataFrame grid_search_rss_varbvsr_sparse(
    const  Eigen::Map<Eigen::SparseMatrix<double> > SiRiS,
    const Rcpp::NumericVector &sigma_beta,
    const Rcpp::NumericVector &logodds,
    const Rcpp::NumericVector &betahat,
    const Rcpp::NumericVector &se,
    const Rcpp::NumericVector &alpha0,
    const Rcpp::NumericVector &mu0,
    const Rcpp::NumericVector &SiRiSr0,
    const double tolerance,
    const int itermax,
    Rcpp::LogicalVector lnz_tol,
    const int n=1,const int grainsize=1){
  
  using namespace Rcpp;
  using namespace tbb;
  using namespace RcppParallel;
  bool islnz_tol = lnz_tol(0);
  
  size_t p=betahat.size();
  size_t tot_size=logodds.size();
  
  
  
  NumericVector sigma_beta_square(tot_size);
  NumericVector sesquare(p);
  NumericVector q(p);
  
  std::transform(sigma_beta.begin(),sigma_beta.end(),sigma_beta_square.begin(),[](double c) {return c*c;});
  std::transform(se.begin(),se.end(),sesquare.begin(),[](double c) {return c*c;});
  std::transform(betahat.begin(),betahat.end(),sesquare.begin(),q.begin(),std::divides<double>());
  
  
  NumericMatrix s=initialize_smat(sesquare,sigma_beta_square);
  NumericMatrix ssrat=initialize_ssratmat(s,sigma_beta_square);
  
  
  //  RMatrix<double> tSiRiS(SiRiS);
  RVector<double> tsigma_beta(sigma_beta);
  RVector<double> tsigma_beta_square(sigma_beta_square);
  RVector<double> tlogodds(logodds);
  RVector<double> tbetahat(betahat);
  RVector<double> tse(se);
  RVector<double> tsesquare(sesquare);
  RVector<double> tq(q);
  
  RMatrix<double> ts(s);
  RMatrix<double> tssrat(ssrat);
  
  
  RVector<double> talpha0(alpha0);
  RVector<double> tmu0(mu0);
  RVector<double> tSiRiSr0(tSiRiSr0);
  
  NumericMatrix alphas(p,tot_size);
  RMatrix<double>talphas(alphas);
  
  NumericMatrix alpha0s(p,tot_size);
  RMatrix<double>talpha0s(alpha0s);
  
  NumericMatrix alpha1s(p,tot_size);
  RMatrix<double>talpha1s(alpha1s);
  
  NumericMatrix mus(p,tot_size);
  RMatrix<double>tmus(mus);
  
  NumericMatrix mu0s(p,tot_size);
  RMatrix<double>tmu0s(mu0s);
  
  NumericMatrix mu1s(p,tot_size);
  RMatrix<double>tmu1s(mu1s);
  
  
  NumericMatrix SiRiSrs(p,tot_size);
  RMatrix<double>tSiRiSrs(SiRiSrs);
  
  
  mdarray ttalpha0(talpha0.begin(),p);
  rassert(ttalpha0.minCoeff()>0 && "alpha0 has all non-negative values");
  
  Fit_wrap talpha(alpha0.begin(),alphas.begin(),alpha0s.begin(),alpha1s.begin(),p,tot_size);
  Fit_wrap tmu(mu0.begin(),mus.begin(),mu0s.begin(),mu1s.begin(),p,tot_size);
  Fit_wrap tSiRiSr(SiRiSr0.begin(),tSiRiSrs.begin(),p,tot_size);
  
  m2darray ttalpha(talpha.fit_p,p,tot_size);
  //  std::cout<<"ttalph.amincoeff: "<<ttalpha.minCoeff()<<std::endl;
  //  std::cout<<"ttalpha0.mincoeff: "<<ttalpha0.minCoeff()<<std::endl;
  rassert(ttalpha.minCoeff()>0 && "alpha has all non-negative values");
  
  
  // NumericVector rlogodds(tot_size);
  // NumericVector rsigma_beta(tot_size);
  NumericVector rfinal_max_err(tot_size);
  IntegerVector rfinal_btnum(tot_size);
  IntegerVector riternum(tot_size);
  NumericVector ralpha_mean(tot_size);
  NumericVector rmu_mean(tot_size);
  NumericVector rpvevec(tot_size);
  NumericVector rlnZ(tot_size);
  
  // std::fill(lnZ.begin(),lnZ.end(),0);
  // std::fill(lnZ0.begin(),lnZ.end(),0);
  
  NumericVector lnZ0(tot_size);
  
  NumericVector lnZ00(tot_size);
  
  NumericVector rvec(tot_size);
  NumericVector mtp(tot_size);
  
  RVector<double> tlnZ0(lnZ0);
  RVector<double> tlnZ00(lnZ00);
  RVector<double> trvec(rvec);
  RVector<double> tmtp(mtp);
  
  
  // RVector<double> trlogodds(rlogodds);
  // RVector<double> trsigma_beta(rsigma_beta);
  RVector<double> trfinal_max_err(rfinal_max_err);
  RVector<int> trfinal_btnum(rfinal_btnum);
  IntegerVector triternum(riternum);
  RVector<double> tralpha_mean(ralpha_mean);
  RVector<double> trmu_mean(rmu_mean);
  RVector<double> trpvevec(rpvevec);
  RVector<double> trlnZ(rlnZ);
  
  
  
  Data_wrap data_wrap(tbetahat.begin(),tse.begin(),sesquare.begin(),q.begin(),SiRiS,p);
  Param_wrap param(tlogodds.begin(),tsigma_beta.begin(),tsigma_beta_square.begin(),ts.begin(),tssrat.begin(),tot_size,p);
  Fit_res fit_res(trlnZ.begin(),
                  tlnZ0.begin(),
                  tlnZ00.begin(),
                  tralpha_mean.begin(),
                  trmu_mean.begin(),
                  trpvevec.begin(),
                  triternum.begin(),
                  trfinal_max_err.begin(),
                  trfinal_btnum.begin(),
                  trvec.begin(),tmtp.begin(),
                  p,
                  tot_size);
  
  std::vector<int> forward_range(p);
  for(size_t i=0;i<p;i++){
    forward_range[i]=i;
  }
  
  std::vector<int> reverse_range(p);
  std::reverse_copy(forward_range.begin(),forward_range.end(),reverse_range.begin());
  
  //  Rcpp::Rcout  <<"Initializing rssr_obj"<<std::endl;
  rssr rssr_obj(&talpha,&tmu,&tSiRiSr,&data_wrap,&param,&fit_res,tolerance,n,itermax,forward_range,reverse_range);
  //Rcpp::Rcout  <<"Starting parallel Algorithm!"<<std::endl;
  //rssr_obj(0,tot_size);
  //Rcpp::Rcout<<"Starting Serial Algorithm!"<<std::endl;
  //  rssr_obj(blocked_range<size_t>(0,tot_size));
  parallel_for(blocked_range<size_t>(0,tot_size,grainsize),rssr_obj);
  //    Rcpp::Rcout<<"Finished!"<<std::endl;
  //  std::cout<<"Checking then Returning!"<<std::endl;
  assert(logodds.size()==sigma_beta.size() &&" logodds==sigma_beta");
  assert(final_max_err.size()==sigma_beta.size() &&" logodds==sigma_beta");
  assert(iternum.size()==sigma_beta.size() &&" logodds==sigma_beta");
  
  
  
  
  
  
  
  assert(std::accumulate(rlogodds.begin(),rlogodds.end(),0.0)!=0 &&"rlogodds shouldn't sum to 0");
  
  assert(std::accumulate(rsigma_beta.begin(),rsigma_beta.end(),0.0)!=0 &&"rsigma_beta shouldn't sum to 0");
  assert(std::accumulate(rfinal_max_err.begin(),rfinal_max_err.end(),0.0)!=0 &&"rfinal_max_err shouldn't sum to 0");
  assert(std::accumulate(riternum.begin(),riternum.end(),0.0)!=0 &&"rfinal_max_err shouldn't sum to 0");
  assert(std::accumulate(ralpha_mean.begin(),ralpha_mean.end(),0.0)!=0 &&"ralpha_mean shouldn't sum to 0");
  assert(std::accumulate(rmu_mean.begin(),rmu_mean.end(),0.0)!=0 &&"rmu_mean shouldn't sum to 0");
  assert(std::accumulate(rmu_mean.begin(),rmu_mean.end(),0.0)!=0 &&"rmu_mean shouldn't sum to 0");
  assert(std::accumulate(rpvevec.begin(),rpvevec.end(),0.0)!=0 &&"rpvevec shouldn't sum to 0");
  assert(std::accumulate(rlnZ.begin(),rlnZ.end(),0.0)!=0 &&"rlnZ shouldn't sum to 0");
  
  //    std::cout<<"Returning!"<<std::endl;
  
  
  return(Rcpp::DataFrame::create(_["logodds"]=logodds,
                                 _["sigb"]=sigma_beta,
                                 _["rel_err"]=rfinal_max_err,
                                 _["iterations"]=riternum,
                                 _["alpha_mean"]=ralpha_mean,
                                 _["mu_mean"]=rmu_mean,
                                 _["pve"]=rpvevec,
                                 _["lnZ"]=rlnZ));
}






