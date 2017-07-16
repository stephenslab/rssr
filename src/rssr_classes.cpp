#include <RcppEigen.h>
#include "rssr.h"

// [[Rcpp::depends(RcppProgress)]]
//[[Rcpp::depends(RcppParallel)]]


template<typename T, typename D> void zero_out(T &dyn_array, const D fillval){
#ifndef EIGEN_NO_DEBUG
  std::fill(dyn_array.local().begin(),dyn_array.local().end(),fillval);
#endif
}
  

using namespace Eigen;

Fit_wrap::Fit_wrap(const double *const ofit_p_,double* fit_,double* fit0_,double* fit1_,const size_t p_, const size_t gridsize_):ofit_p(ofit_p_),
																 p(p_),
																 gridsize(gridsize_),
																 fit_p(fit_),
																 fit0_p(fit0_),
																 fit1_p(fit1_){
  
  
  c_mdarray ofit(ofit_p,p);
  m2darray fit(fit_p,p,gridsize);
  for(size_t i=0; i<gridsize;i++){
    fit.col(i)=ofit;
  }
  
  if(ofit.minCoeff()>0){
    //    std::cout<<"oarray.minCoeff: "<<ofit.minCoeff()<<std::endl;
    //    std::cout<<"fit.mincoeff: "<<fit.minCoeff()<<std::endl;
    rassert(fit.minCoeff()>0 && "if ofit_p is greater than zero, so should the fit matrix");
  }
  //  fit.colwise()=ofit;
  
  //  rassert((fit.col(0)==ofit) &&"fit col 0 should equal alpha0");
#ifndef NDEBUG
  rassert(ofit.size()==p &&"fit should have size of number of SNPs (p)");
  rassert(ofit.size()==p*gridsize &&"fit should have size of number of SNPs (p)");
  rassert(isnan(fit).sum()==0 && "fit has no NaN  values");
#endif
																 }

Fit_wrap::Fit_wrap(const double* const ofit_p_,double* fit_,const size_t p_,const size_t gridsize_):ofit_p(ofit_p_),
												   p(p_),
												  gridsize(gridsize_),
												fit_p(fit_),
												fit0_p(NULL),
												fit1_p(NULL){
    c_mdarray ofit(ofit_p,p);
    m2darray fit(fit_p,p,gridsize);
    for(size_t i=0; i<gridsize;i++){
      fit.col(i)=ofit;
    }
    //    fit.colwise()=ofit;
    
    
#ifndef NDEBUG
    rassert(ofit.size()==p &&"fit should have size of number of SNPs (p)");
    rassert(ofit.size()==p*gridsize &&"fit should have size of number of SNPs (p)");
    rassert(isnan(fit).sum()==0 && "fit has no NaN  values");
#endif
  }														       

Data_wrap::Data_wrap(const double* betahat_,const double* se_,const double* sesquare_,const double* q_,
	  const sparseMatrix_external siris_, const size_t p_):p(p_),
							       isSparse(true),
								 betahat_p(betahat_),
								 se_p(se_),
								 siris_p(NULL),
								 msiris(siris_),
								 sesquare_p(sesquare_),
								 q_p(q_){
  c_mdarray betahat (betahat_p,p);
  c_mdarray se (se_p,p);
  //  c_m2darray siris (siris_p,p,p);
  c_mdarray sesquare (sesquare_p,p);
  c_mdarray q(q_p,p);
  
  
  rassert(sesquare.cols()==1 && "sesquare should have 1 column");
  rassert(se.rows()==sesquare.rows() && "se should be same length as sesquare");
  rassert(se.size()==p && "se should have p elements");
  rassert(q.size()==p && "q should have p elements");
  rassert(betahat.size()==p && "q should have p elements");
  
  rassert(se.minCoeff()>0 && "se has all non-negative values");
  rassert(sesquare.minCoeff()>0 && "sesquare has all non-negative values");
  rassert(isnan(sesquare).sum()==0 && "sesquare has no NaN  values");
  rassert(isnan(q).sum()==0 && "q has no NaN  values");
  rassert(isnan(betahat).sum()==0 && "betahat has no NaN  values");


    
								 }


Data_wrap::Data_wrap(const double* betahat_,const double* se_,const double* sesquare_,const double* q_,const double *siris_, const size_t p_):p(p_),
																	      isSparse(false),
																	      betahat_p(betahat_),
																	      se_p(se_),
																	      siris_p(siris_),
  msiris(sparseMatrix_external(p_,p_,0,NULL,NULL,NULL)),												      sesquare_p(sesquare_),
																	      q_p(q_){
  c_mdarray betahat (betahat_p,p);
  c_mdarray se (se_p,p);
  c_m2darray siris (siris_p,p,p);
  c_mdarray sesquare (sesquare_p,p);
  c_mdarray q(q_p,p);
  
  
  rassert(sesquare.cols()==1 && "sesquare should have 1 column");
  rassert(se.rows()==sesquare.rows() && "se should be same length as sesquare");
  rassert(se.size()==p && "se should have p elements");
  rassert(q.size()==p && "q should have p elements");
  rassert(betahat.size()==p && "q should have p elements");
  
  rassert(se.minCoeff()>0 && "se has all non-negative values");
  rassert(sesquare.minCoeff()>0 && "sesquare has all non-negative values");
  rassert(isnan(sesquare).sum()==0 && "sesquare has no NaN  values");
  rassert(isnan(q).sum()==0 && "q has no NaN  values");
  rassert(isnan(betahat).sum()==0 && "betahat has no NaN  values");


    
																	      }

Param_wrap::Param_wrap(const double* logodds_,const double* sigb_,const double* sigma_beta_square_,const double* s_, const double* ssrat_,const size_t gridsize_,const size_t p_):gridsize(gridsize_),
																	   p(p_),
																	   logodds_p(logodds_),
																	   sigma_beta_p(sigb_),
																	   sigma_beta_square_p(sigma_beta_square_),
																	   s_p(s_),
																	   ssrat_p(ssrat_){
  

  c_mdarray logodds(logodds_p,gridsize);
  c_mdarray sigma_beta(sigma_beta_p,gridsize);
  c_mdarray sigma_beta_square(sigma_beta_square_p,gridsize);
  c_m2darray s(s_p,p,gridsize);
  c_m2darray ssrat(ssrat_p,p,gridsize);
  
  
  
  
  rassert(sigma_beta.size()==sigma_beta_square.size() && "sigma_beta should be same length as sigma_beta_square");
  rassert(logodds.size()==gridsize && "logodds should have size of gridsize");
  
  rassert(s.rows()==p && "s should have p rows");
  rassert(ssrat.rows()==p && "ssrat should have p rows");

  rassert(s.cols()==gridsize && "s should have gridsize cols");
  rassert(ssrat.cols()==gridsize && "ssrat should have gridsize cols");

  rassert(isnan(s).sum()==0 && "s has no NaN  values");
  rassert(isnan(ssrat).sum()==0 && "ssrat has no NaN  values");
  
  rassert(isnan(sigma_beta).sum()==0 && "sigma_beta has no NaN  values");
  rassert(isnan(sigma_beta_square).sum()==0 && "sigma_beta_square has no NaN  values");
  // rassert(isnan(ssrat).sum()==0 && "ssrat has no NaN  values");
  rassert(isnan(logodds).sum()==0 && "logodds has no NaN  values");
}


Param_wrap::Param_wrap(const double* sigb_,const double* sigma_beta_square_,const double* s_, const double* ssrat_,const size_t gridsize_,const size_t p_):gridsize(gridsize_),
																			   p(p_),
																			   logodds_p(NULL),
																			   sigma_beta_p(sigb_),
																			   sigma_beta_square_p(sigma_beta_square_),
																			   s_p(s_),
																			   ssrat_p(ssrat_){
  
  
  c_mdarray sigma_beta(sigma_beta_p,gridsize);
  c_mdarray sigma_beta_square(sigma_beta_square_p,gridsize);
  c_m2darray s(s_p,p,gridsize);
  c_m2darray ssrat(ssrat_p,p,gridsize);
  
  
  rassert(sigma_beta.size()==sigma_beta_square.size() && "sigma_beta should be same length as sigma_beta_square");
  rassert(s.rows()==p && "s should have p rows");
  rassert(ssrat.rows()==p && "ssrat should have p rows");
  rassert(s.cols()==gridsize && "s should have gridsize cols");
  rassert(ssrat.cols()==gridsize && "ssrat should have gridsize cols");
  rassert(isnan(s).sum()==0 && "s has no NaN  values");
  rassert(isnan(ssrat).sum()==0 && "ssrat has no NaN  values");
  rassert(isnan(sigma_beta).sum()==0 && "sigma_beta has no NaN  values");
  rassert(isnan(sigma_beta_square).sum()==0 && "sigma_beta_square has no NaN  values");
																	   }








  

Fit_res::Fit_res(double *const lnZ_,
		 double *const lnZ0_,
		 double *const lnZ00_,
		 double *const alpha_mean_,
		 double *const mu_mean_,
		 double *const pvevec_,
		 int *const  iternum_,
		 double *const final_max_err_,
		 int *const final_btnum_,
		 double *const rvec_,
		 double *const mtp_,const size_t p_,const size_t gridsize_):  lnZ_p(lnZ_),
									      lnZ0_p(lnZ0_),
									      lnZ00_p(lnZ00_),
									      alpha_mean_p(alpha_mean_),
									      mu_mean_p(mu_mean_),
									      pvevec_p(pvevec_),
									      iternum_p(iternum_),
									      final_max_err_p(final_max_err_),
									      final_btnum_p(final_btnum_),
									      rvec_p(rvec_),
									      mtp_p(mtp_),
									      p(p_),
									      gridsize(gridsize_)
{
  


  
  
  
}

    
    
    
  
  

rssr::rssr(Fit_wrap *alphas_,
	   Fit_wrap *mus_,
	   Fit_wrap *sirisrs_,
	   const Data_wrap *datas_,
	   const Param_wrap *params_,
	   Fit_res *results_,
	   const double tolerance_,
	   const size_t n_,
	   const size_t itermax_,
	   const std::vector<int> &forward_range_,
	   const std::vector<int> &reverse_range_): isSparse(datas_->isSparse),
						    gridsize(params_->gridsize),
						    p(mus_->p),
						    alpha(alphas_->fit_p,p,gridsize),
						    alpha0(alphas_->fit0_p,p,gridsize),
						    alpha1(alphas_->fit1_p,p,gridsize),
						    mu(mus_->fit_p,p,gridsize),
						    mu0(mus_->fit0_p,p,gridsize),
						    mu1(mus_->fit1_p,p,gridsize),
						    SiRiSr(sirisrs_->fit_p,p,gridsize),
						    n(n_),
						    q(datas_->q_p,p),
						    logodds(params_->logodds_p,gridsize),
						    sigma_beta(params_->sigma_beta_p,gridsize),
						    sigma_beta_square(params_->sigma_beta_square_p,gridsize),
						    s(params_->s_p,p,gridsize),
						    ssrat(params_->ssrat_p,p,gridsize),
						    sesquare(datas_->sesquare_p,p),
						    itermax(itermax_),
						    tolerance(tolerance_),
						    forward_range(forward_range_),
						    reverse_range(reverse_range_),
						    betahat(datas_->betahat_p,p),
						    se(datas_->se_p,p),
						    SiRiS(datas_->siris_p,p,p),
						    sSiRiS(datas_->msiris),
						    lnZ(results_->lnZ_p,gridsize),
						    lnZ0(results_->lnZ0_p,gridsize),
						    lnZ00(results_->lnZ00_p,gridsize),
						    alpha_mean(results_->alpha_mean_p,gridsize),
						    mu_mean(results_->mu_mean_p,gridsize),
						    pvevec(results_->pvevec_p,gridsize),
						    iternum(results_->iternum_p,gridsize),
						    final_btnum(results_->final_btnum_p,gridsize),
						    final_max_err(results_->final_max_err_p,gridsize),
						    rvec(results_->rvec_p,gridsize),
						    mtp(results_->mtp_p,gridsize)
{
  
  rassert(alpha.minCoeff()>0 && "alpha has all non-negative values");
  rassert(isnan(s).sum()==0 && "s has no NaN  values");
  rassert(isnan(sigma_beta).sum()==0 && "sigma_beta has no NaN  values");
  rassert(isnan(sigma_beta_square).sum()==0 && "sigma_beta_square has no NaN  values");
  rassert(isnan(ssrat).sum()==0 && "ssrat has no NaN  values");
  rassert(isnan(logodds).sum()==0 && "logodds has no NaN  values");
  
  
  rassert(isnan(alpha).sum()==0 && "alpha has no NaN  values");
  rassert(isnan(mu).sum()==0 && "mu has no NaN  values");
  rassert(isnan(SiRiSr).sum()==0 && "SiRiSr has no NaN  values");
  rassert(isnan(alpha0).sum()==0 && "alpha0 has no NaN values");
  rassert(isnan(alpha1).sum()==0 && "alpha1 has no NaN values");
  rassert(isnan(mu0).sum()==0 && "mu0 has no NaN  values");
  rassert(isnan(mu1).sum()==0 && "mu1 has no NaN  values");
  

  rassert(isnan(alpha).sum()==0 && "alpha has no NaN starting values");
  rassert(isnan(mu).sum()==0 && "mu has no NaN starting values");
  rassert(isnan(SiRiSr).sum()==0 && "SiRiSr has no NaN starting values"); 
}



void rssr::calculate_mtp(const blocked_range <size_t > &r)const {

  for(size_t i=r.begin();i!=r.end();i++){
    mtp(i)=
      -std::sqrt(((alpha1.col(i)-alpha0.col(i)).square()).sum()+((mu1.col(i)-mu0.col(i)).square()).sum())/std::sqrt(((alpha.col(i)-2*alpha1+alpha0).square()).sum()+((mu.col(i)-2*mu1.col(i)+mu0.col(i)).square()).sum()+double_lim::epsilon());
  }
}


void rssr::calc_lnZ(const blocked_range<size_t> & r)const {
  size_t rsize=r.end()-r.begin();

  //  std::cout<<"(lnZ) rsize is :"<<rsize<<std::endl;
  rassert(lnZ.size()>=rsize && "lnZ.size() is smaller than chunk ");


  ////std::cout<<"First check that alpha is ready to go"<<std::endl;
  rassert(alpha.cols()>=rsize && "alpha has to have colno>=rsize");
  rassert(alpha.rows()==p && "alpha has to have p rows");
  rassert(alpha(0,0)!=2 && "can access alpha");
  rassert(isnan(alpha).sum()==0 && "alpha is all non NaN");

  ////std::cout<<"Now check mu"<<std::endl;
  rassert(mu.cols()>=rsize && "mu has to have colno>=rsize");
  rassert(mu.rows()==p && "mu has to have p rows");
  rassert(mu(0,0)!=2 && "can access mu");
  rassert(isnan(mu).sum()==0 && "mu is all non NaN");
  
  ////std::cout<<"Now check SiRiSr"<<std::endl;
  rassert(SiRiSr.cols()>=rsize && "SiRiSr has to have colno>=rsize");
  rassert(SiRiSr.rows()==p && "SiRiSr has to have p rows");
  rassert(SiRiSr(0,0)!=2 && "can access SiRiSr");
  rassert(isnan(SiRiSr).sum()==0 && "SiRiSr is all non NaN");

  ////std::cout<<"Now check s"<<std::endl;
  rassert(s.cols()>=rsize && "s has to have colno>=rsize");
  rassert(s.rows()==p && "s has to have p rows");

  rassert(s.sum()!=0 && "can access elements of  s");
  rassert(isnan(s).sum()==0 && "s is all non NaN");


  rassert(isnan(q).sum()==0 && "q is all non NaN");
  ////std::cout<<"Now check ssrat"<<std::endl;
  rassert(ssrat.cols()>=rsize && "ssrat has to have colno>=rsize");
  rassert(ssrat.rows()==p && "ssrat has to have p rows");
  rassert(ssrat(0,0)!=2 && "can access ssrat");
  rassert(isnan(ssrat).sum()==0 && "ssrat is all non NaN");

  ////std::cout<<"Now check sesquare"<<std::endl;
  rassert(sesquare.size()==p && "sesquare has to have p elem");
  rassert(sesquare(0)!=2 && "can access sesquare");
  rassert(isnan(sesquare).sum()==0 && "sesquare is all non NaN");

  ////std::cout<<"Now check logodds"<<std::endl;
  rassert(logodds.size()>=rsize && "logodds has to have rsize elem");
  rassert(logodds(0)!=2 && "can access logodds");
  rassert(isnan(logodds).sum()==0 && "logodds is all non NaN");

  ////std::cout<<"Now check sigma_beta"<<std::endl;
  rassert(sigma_beta.size()>=rsize && "sigma_beta has to have rsize elem");
  rassert(sigma_beta(0)!=2 && "can access sigma_beta");
  rassert(isnan(sigma_beta).sum()==0 && "sigma_beta is all non NaN");
  
  ////std::cout<<"about to start to calculate lnZ"<<std::endl;
  for(size_t j=r.begin(); j!=r.end();j++){

    lnZ(j)=calculate_lnZ(q,
			 alpha.col(j)*mu.col(j),
			 SiRiSr.col(j),
			 logodds(j),
			 sesquare,
			 alpha.col(j),
			 mu.col(j),
			 s.col(j),
			 sigma_beta(j));
    if(isnan(lnZ(j))){
	std::cerr<<j<<std::endl;
	std::cerr<<"lnZ(j):"<<lnZ(j)<<std::endl;
	std::cerr<<"sigma_beta(j):"<<sigma_beta(j)<<std::endl;
	std::cerr<<"logodds(j):"<<logodds(j)<<std::endl;
	
	rassert(!isnan(lnZ(j))&& "lnZ(j) is finite ");
      }
  }
}

  

void rssr::rss_vb_iter(const blocked_range<size_t> r,bool reverse)const{
  if(isSparse){
    rss_vb_iter_sparse(r,reverse);
  }else{ 
    size_t rstart=r.begin();
    size_t rsize= r.end()-rstart;			  
    if(reverse){
      for(const int& i: reverse_range){
	for(size_t j=r.begin();j!=r.end();j++){
	  rvec(j)=alpha(i,j)*mu(i,j);
	  mu(i,j)=s(i,j)*(q(i)+rvec(j)/sesquare(i)-SiRiSr(i,j));
	  alpha(i,j)=sigmoid(logodds(j)+0.5*(ssrat(i,j)+(mu(i,j)*mu(i,j))/s(i,j)));
	  SiRiSr.col(j)+=SiRiS.col(i).array()*(alpha(i,j)*mu(i,j)-rvec(j));
	}
      }
    }else{
      for(const int& i: forward_range){
	for(size_t j=r.begin();j!=r.end();j++){
	  rvec(j)=alpha(i,j)*mu(i,j);
	  mu(i,j)=s(i,j)*(q(i)+rvec(j)/sesquare(i)-SiRiSr(i,j));
	  alpha(i,j)=sigmoid(logodds(j)+0.5*(ssrat(i,j)+(mu(i,j)*mu(i,j))/s(i,j)));
	  SiRiSr.col(j)+=SiRiS.col(i).array()*(alpha(i,j)*mu(i,j)-rvec(j));
	}	    
      }
    }
    rassert(isnan(alpha).sum()==0 && "alpha has no NaN  values");
    rassert(isnan(mu).sum()==0 && "mu has no NaN  values");
    rassert(isnan(SiRiSr).sum()==0 && "SiRiSr has no NaN  values");
  }
}

void rssr::rss_vb_iter_sparse(const blocked_range<size_t> r,bool reverse)const{
    
  size_t rstart=r.begin();
  size_t rsize= r.end()-rstart;
  const std::vector<int>* grange;
  if(reverse){
    grange=&reverse_range;
  }else{
    grange=&forward_range;
  }
  for(const int& i: *grange){
    for(size_t j=r.begin();j!=r.end();j++){
      rvec(j)=mu(i,j);
      mu(i,j)=s(i,j)*(q(i)+rvec(j)/sesquare(i)-SiRiSr(i,j));
      alpha(i,j)=sigmoid(logodds(j)+0.5*(ssrat(i,j)+(mu(i,j)*mu(i,j))/s(i,j)));
      SiRiSr.col(j)+=sSiRiS.col(i)*(alpha(i,j)*mu(i,j)-rvec(j));
    }	    
  }
  rassert(isnan(mu).sum()==0 && "mu has no NaN  values");
  rassert(isnan(SiRiSr).sum()==0 && "SiRiSr has no NaN  values");	
}
  
    
  
void rssr::squarem_step(const blocked_range<size_t> & r) const{
  if(isSparse){
    squarem_step_sparse(r);
  }else{
    calculate_mtp(r);
    for(size_t j=r.begin(); j!=r.end();j++){
      double tmtp=mtp(j);
      if(tmtp>(-1)){
	tmtp=-1;
	mtp(j) =-1;
      }else{
	if(!std::isnan(tmtp)){
	  //<<"alpha_step:"<<std::endl;
	  alpha.col(j)=alpha0.col(j)-2*tmtp*(alpha1.col(j)-alpha0.col(j))+(tmtp*tmtp)*((alpha.col(j)-2*alpha1.col(j)+alpha0.col(j)));
	  //<<"mu_step:"<<std::endl;
	  mu.col(j)=mu0.col(j)-2*tmtp*((mu1.col(j)-mu0.col(j)))+(tmtp*tmtp)*(mu.col(j)-2*mu1.col(j)+mu0.col(j));
	  //<<"SiRiSr_step:"<<std::endl;
	  SiRiSr.col(j)=(SiRiS.matrix()*(alpha.col(j)*mu.col(j)).matrix()).array();
	}
      }
    }
    rassert(isnan(alpha).sum()==0 && "alpha has no NaN  values");
    rassert(isnan(mu).sum()==0 && "mu has no NaN  values");
    rassert(isnan(SiRiSr).sum()==0 && "SiRiSr has no NaN  values");
  }
}

void rssr::squarem_step_sparse(const blocked_range<size_t> & r) const{    
  calculate_mtp(r);
  for(size_t j=r.begin(); j!=r.end();j++){
    double tmtp=mtp(j);
    if(tmtp>(-1)){
      tmtp=-1;
      mtp(j) =-1;
    }else{
      if(!std::isnan(tmtp)){
	//<<"alpha_step:"<<std::endl;
	alpha.col(j)=alpha0.col(j)-2*tmtp*(alpha1.col(j)-alpha0.col(j))+(tmtp*tmtp)*((alpha.col(j)-2*alpha1.col(j)+alpha0.col(j)));
	//<<"mu_step:"<<std::endl;
	mu.col(j)=mu0.col(j)-2*tmtp*((mu1.col(j)-mu0.col(j)))+(tmtp*tmtp)*(mu.col(j)-2*mu1.col(j)+mu0.col(j));
	//<<"SiRiSr_step:"<<std::endl;
	SiRiSr.col(j)=(sSiRiS*(alpha.col(j)*mu.col(j)).matrix());
      }
    }
  }
  rassert(isnan(alpha).sum()==0 && "alpha has no NaN  values");
  rassert(isnan(mu).sum()==0 && "mu has no NaN  values");
  rassert(isnan(SiRiSr).sum()==0 && "SiRiSr has no NaN  values");
}


      
void rssr::squarem_backtrack(const blocked_range<size_t> & r)const {
  if(isSparse){
    squarem_backtrack_sparse(r);
  }else{
    for(size_t j=r.begin();j!=r.end();j++){
      double tmtp=mtp(j);
      double tlnZ=lnZ(j);
      double tlnZ0=lnZ0(j);
      if(tmtp<(-1) && (tlnZ < tlnZ0)){
	size_t num_bt=0;
	while((tlnZ<tlnZ0) && (num_bt < max_backtrack)){
	  tmtp=0.5*(tmtp-1);
	  alpha.col(j) = alpha0.col(j)-2*tmtp*(alpha1.col(j)-alpha0.col(j))+(tmtp*tmtp)*(alpha.col(j)-2*alpha1.col(j)+alpha0.col(j));
	  mu.col(j) = mu0.col(j)-2*tmtp*(mu1.col(j)-mu0.col(j))+(tmtp*tmtp)*(mu.col(j)-2*mu1.col(j)+mu0.col(j));
	  SiRiSr.col(j) = (SiRiS.matrix()*(alpha.col(j)*mu.col(j)).matrix()).array();
	  rss_vb_iter(blocked_range<size_t>(j,j),false);
	  rss_vb_iter(blocked_range<size_t>(j,j),true);
	      
	  calc_lnZ(blocked_range<size_t>(j,j));
	  num_bt=num_bt+1;
	  final_btnum(j)+=1;
	}
	if(num_bt==max_backtrack){
	  alpha.col(j)=alpha0.col(j);
	  mu.col(j)=mu0.col(j);
	  SiRiSr.col(j) = (SiRiS.matrix()*(alpha.col(j)*mu.col(j)).matrix()).array();
	  calc_lnZ(blocked_range<size_t>(j,j));
	  lnZ0(j)=lnZ00(j);
	}
      }
      mtp(j)=tmtp;
    }
  }
}


void rssr::squarem_backtrack_sparse(const blocked_range<size_t> & r)const {

  for(size_t j=r.begin();j!=r.end();j++){
    double tmtp=mtp(j);
    double tlnZ=lnZ(j);
    double tlnZ0=lnZ0(j);
    if(tmtp<(-1) && (tlnZ < tlnZ0)){
      size_t num_bt=0;
      while((tlnZ<tlnZ0) && (num_bt < max_backtrack)){
	tmtp=0.5*(tmtp-1);
	alpha.col(j) = alpha0.col(j)-2*tmtp*(alpha1.col(j)-alpha0.col(j))+(tmtp*tmtp)*(alpha.col(j)-2*alpha1.col(j)+alpha0.col(j));
	mu.col(j) = mu0.col(j)-2*tmtp*(mu1.col(j)-mu0.col(j))+(tmtp*tmtp)*(mu.col(j)-2*mu1.col(j)+mu0.col(j));
	SiRiSr.col(j) = (sSiRiS*(alpha.col(j)*mu.col(j)).matrix());
	rss_vb_iter(blocked_range<size_t>(j,j),false);
	rss_vb_iter(blocked_range<size_t>(j,j),true);
	      
	calc_lnZ(blocked_range<size_t>(j,j));
	num_bt=num_bt+1;
	final_btnum(j)+=1;
      }
      if(num_bt==max_backtrack){
	alpha.col(j)=alpha0.col(j);
	mu.col(j)=mu0.col(j);
	SiRiSr.col(j) = (sSiRiS*(alpha.col(j)*mu.col(j)).matrix());
	calc_lnZ(blocked_range<size_t>(j,j));
	lnZ0(j)=lnZ00(j);
      }
    }
    mtp(j)=tmtp;
  }
}
    
void rssr::calc_err(const blocked_range<size_t> &r,double &max_err)const {
      
  max_err=0;
  if(use_lnztol){
	
    for(size_t j=r.begin(); j!=r.end();j++){
      double tmax_err=rel_err(lnZ(j),lnZ0(j));
      if(tmax_err>max_err){
	max_err=tmax_err;
      }
    }
  }else{
    for(size_t j=r.begin(); j!=r.end();j++){
      double tmax_err=find_maxerr(alpha.col(j),alpha0.col(j),alpha.col(j)*mu.col(j),alpha0.col(j)*mu0.col(j));
      if(tmax_err>max_err){
	max_err=tmax_err;
      }
    }
  }
}
void rssr::compute_pve(const blocked_range<size_t> &r)const {
  // size_t rstart=r.begin();
  // size_t rsize= r.end()-rstart;
  for(size_t j=r.begin();j!=r.end();j++){
    pvevec(j)=(SiRiSr.col(j).matrix().transpose()*(alpha.col(j)*mu.col(j)).matrix());
    pvevec(j)/=(double)n;
  }
}    

  
void rssr::operator()(    const blocked_range<size_t> &r)const {

    
  double max_err=1;
  size_t rstart=r.begin();
  size_t rsize= r.end()-rstart;
  //  blocked_range<size_t> r(0,rsize);
  //  tls_init(r);
  calc_lnZ(r);
  int iter=0;
  
  // Rcpp::Rcout<<"Initializing lnZ Segment "<<std::endl;    
  calc_lnZ(r);
  //Rcpp::Rcout<<"Initializing lnZ0 Segment "<<std::endl;
  lnZ0.segment(rstart,rsize)=lnZ.segment(rstart,rsize);
  while(max_err>tolerance){
    
    //Rcpp::Rcout<<"Setting lnZ00 Segment "<<std::endl;
    lnZ00.segment(rstart,rsize)=lnZ0.segment(rstart,rsize);
    //Rcpp::Rcout<<"Setting lnZ00 Segment "<<std::endl;
    lnZ0.segment(rstart,rsize)=lnZ.segment(rstart,rsize);
    //Rcpp::Rcout<<"Setting alpha0 block "<<std::endl;
    alpha0.block(0,rstart,p,rsize)=alpha.block(0,rstart,p,rsize);
    //Rcpp::Rcout<<"Setting mu0 block "<<std::endl;
    mu0.block(0,rstart,p,rsize)=mu.block(0,rstart,p,rsize);
    
    bool reverse = iter%2!=0;
    //Rcpp::Rcout<<"Performing first Iteration "<<std::endl;
    rss_vb_iter(r,reverse);
    //Rcpp::Rcout<<"Setting alpha1 block "<<std::endl;
    alpha1.block(0,rstart,p,rsize)=alpha.block(0,rstart,p,rsize);
    mu1.block(0,rstart,p,rsize)=mu.block(0,rstart,p,rsize);
    //Rcpp::Rcout<<"Performing second Iteration "<<std::endl;
    rss_vb_iter(r,reverse);
    
    rassert(isnan(alpha.block(0,rstart,p,rsize)).sum()==0 && "alpha has no NaN  values");
    rassert(isnan(mu.block(0,rstart,p,rsize)).sum()==0 && "mu has no NaN  values");
    rassert(isnan(SiRiSr.block(0,rstart,p,rsize)).sum()==0 && "SiRiSr has no NaN  values");
    
    rassert(isnan(alpha0.block(0,rstart,p,rsize)).sum()==0 && "alpha0 has no NaN values");
    rassert(isnan(alpha1.block(0,rstart,p,rsize)).sum()==0 && "alpha1 has no NaN values");
    


    rassert(isnan(mu0.block(0,rstart,p,rsize)).sum()==0 && "mu0 has no NaN  values");
    rassert(isnan(mu1.block(0,rstart,p,rsize)).sum()==0 && "mu1 has no NaN  values");

    
    rassert(isnan(lnZ.segment(rstart,rsize)).sum()==0 && "lnZ has no NaN  values");
    rassert(isnan(lnZ0.segment(rstart,rsize)).sum()==0 && "lnZ0 has no NaN  values");

    //Rcpp::Rcout<<"Squarem Step "<<std::endl;
    squarem_step(r);
    //Rcpp::Rcout<<"Third Iteration "<<std::endl;
    rss_vb_iter(r,reverse);
    calc_lnZ(r);
    //Rcpp::Rcout<<"Squarem Backtrack "<<std::endl;
    squarem_backtrack(r);
    //Rcpp::Rcout<<"Error Calculation "<<std::endl;
    calc_err(r,max_err);
    if(iter>itermax){
      for(size_t  j=r.begin();j!=r.end();j++){
	final_max_err(j)=max_err;
	iternum(j)=iter;
      }
      compute_pve(r);
      //      results->write_res(qr);
      break;
      
    }
    iter++;
  }
  //Rcpp::Rcout<<"Error Calculation "<<std::endl;
    
  for(size_t  j=r.begin();j!=r.end();j++){
    final_max_err(j)=max_err;
    iternum(j)=iter;
    alpha_mean(j)=alpha.col(j).mean();
    mu_mean(j)=mu.col(j).mean(); 
  }
  
  rassert(isnan(alpha_mean).sum()==0 && "alpha_mean is entirely finite");
  //  rassert(1==2 && "just checking that rasserts actually do something" && //std::cout<<"Something happened"<<std::endl);

  
  compute_pve(r);
  //  results->write_res(qr);
}
    
  
rssr_norm::rssr_norm(Fit_wrap *mus_,
		     Fit_wrap *sirisrs_,
		     const Data_wrap *datas_,
		     Param_wrap *params_,
		     Fit_res *results_,
		     const double tolerance_, const size_t n_,
		     const size_t itermax_,
		     const std::vector<int> &forward_range_,
		     const std::vector<int> &reverse_range_
		     ):isSparse(datas_->isSparse),
		       p(mus_->p),
		       gridsize(params_->gridsize),
		       mu(mus_->fit_p,p,gridsize),
		       mu0(mus_->fit0_p,p,gridsize),
		       mu1(mus_->fit1_p,p,gridsize),
		       SiRiSr(sirisrs_->fit_p,p,gridsize),
		       n(n_),
		       q(datas_->q_p,p),
		       sigma_beta(params_->sigma_beta_p,gridsize),
		       sigma_beta_square(params_->sigma_beta_square_p,gridsize),
		       s(params_->s_p,p,gridsize),
		       ssrat(params_->ssrat_p,p,gridsize),
		       sesquare(datas_->sesquare_p,p),
		       itermax(itermax_),
		       tolerance(tolerance_),
		       forward_range(forward_range_),
		       reverse_range(reverse_range_),
		       betahat(datas_->betahat_p,p),
		       se(datas_->se_p,p),
		       SiRiS(datas_->siris_p,p,p),
		       sSiRiS(datas_->msiris),
		       lnZ(results_->lnZ_p,gridsize),
		       lnZ0(results_->lnZ0_p,gridsize),
		       lnZ00(results_->lnZ00_p,gridsize),
		       alpha_mean(results_->alpha_mean_p,gridsize),
		       mu_mean(results_->mu_mean_p,gridsize),
		       pvevec(results_->pvevec_p,gridsize),
		       iternum(results_->iternum_p,gridsize),
		       final_btnum(results_->final_btnum_p,gridsize),
		       final_max_err(results_->final_max_err_p,gridsize),
		       rvec(results_->rvec_p,gridsize),
		       mtp(results_->mtp_p,gridsize){

  rassert(isnan(s).sum()==0 && "s has no NaN  values");
  rassert(isnan(sigma_beta).sum()==0 && "sigma_beta has no NaN  values");
  rassert(isnan(sigma_beta_square).sum()==0 && "sigma_beta_square has no NaN  values");
  rassert(isnan(ssrat).sum()==0 && "ssrat has no NaN  values");
  rassert(isnan(mu).sum()==0 && "mu has no NaN  values");
  rassert(isnan(SiRiSr).sum()==0 && "SiRiSr has no NaN  values");
  
  rassert(isnan(mu0).sum()==0 && "mu0 has no NaN  values");
  rassert(isnan(mu1).sum()==0 && "mu1 has no NaN  values");
  

  rassert(isnan(mu).sum()==0 && "mu has no NaN starting values");
  rassert(isnan(SiRiSr).sum()==0 && "SiRiSr has no NaN starting values"); 

}





void rssr_norm::calculate_mtp(const blocked_range <size_t > &r)const {

  for(size_t i=r.begin();i!=r.end();i++){
    mtp(i)=-std::sqrt(((mu1.col(i)-mu0.col(i)).square()).sum())/std::sqrt(((mu.col(i)-2*mu1.col(i)+mu0.col(i)).square()).sum()+double_lim::epsilon());
  }
}


void rssr_norm::calc_lnZ(const blocked_range<size_t> & r)const {
  size_t rsize=r.end()-r.begin();

    
  rassert(lnZ.size()>=rsize && "lnZ.size() is smaller than chunk ");
  rassert(mu.cols()>=rsize && "mu has to have colno>=rsize");
  rassert(SiRiSr.cols()>=rsize && "SiRiSr has to have colno>=rsize");
  rassert(s.cols()>=rsize && "s has to have colno>=rsize");
  rassert(sigma_beta.size()>=rsize && "sigma_beta has to have size>=rsize");
  rassert(sesquare.size()==p && "sesquare has to have p elements");

  for(size_t j=r.begin(); j!=r.end();j++){
      lnZ(j)=calculate_lnZ(q,
				 mu.col(j),
				 SiRiSr.col(j),
				 sesquare,
				 mu.col(j),
				 s.col(j),
				 sigma_beta(j));
  }
}

  

void rssr_norm::rss_vb_iter(const blocked_range<size_t> r,bool reverse)const{

  if(isSparse){
    rss_vb_iter_sparse(r,reverse);
  }else{
    size_t rstart=r.begin();
    size_t rsize= r.end()-rstart;

    // for(size_t j=r.begin(); j!=r.end();j++){
    //        bool reverse = iternum[j]%2!=0;				  
    if(reverse){
      for(const int& i: reverse_range){
	for(size_t j=r.begin();j!=r.end();j++){
	  rvec(j)=mu(i,j);
	  mu(i,j)=s(i,j)*(q(i)+rvec(j)/sesquare(i)-SiRiSr(i,j));
	  SiRiSr.col(j)+=SiRiS.col(i).array()*(mu(i,j)-rvec(j));
	}
      }
      
    }else{
      for(const int& i: forward_range){
	for(size_t j=r.begin();j!=r.end();j++){
	  rvec(j)=mu(i,j);
	  mu(i,j)=s(i,j)*(q(i)+rvec(j)/sesquare(i)-SiRiSr(i,j));
	  SiRiSr.col(j)+=SiRiS.col(i).array()*(mu(i,j)-rvec(j));
	}	    
      }
    }
    
    rassert(isnan(mu).sum()==0 && "mu has no NaN  values");
    rassert(isnan(SiRiSr).sum()==0 && "SiRiSr has no NaN  values");
    
  }
}
void rssr_norm::rss_vb_iter_sparse(const blocked_range<size_t> r,bool reverse)const{
    
  size_t rstart=r.begin();
  size_t rsize= r.end()-rstart;
  const std::vector<int>* grange;
  if(reverse){
    grange=&reverse_range;
  }else{
    grange=&forward_range;
  }
    

  // for(size_t j=r.begin(); j!=r.end();j++){
  //        bool reverse = iternum[j]%2!=0;
  for(const int& i: *grange){
    for(size_t j=r.begin();j!=r.end();j++){
      rvec(j)=mu(i,j);
      mu(i,j)=s(i,j)*(q(i)+rvec(j)/sesquare(i)-SiRiSr(i,j));
      SiRiSr.col(j)+=sSiRiS.col(i)*(mu(i,j)-rvec(j));
    }	    
  }
  rassert(isnan(mu).sum()==0 && "mu has no NaN  values");
  rassert(isnan(SiRiSr).sum()==0 && "SiRiSr has no NaN  values");
	
}



    
    
  
void rssr_norm::squarem_step(const blocked_range<size_t> & r) const{
  if(isSparse){
    squarem_step_sparse(r);
  }else{    
    calculate_mtp(r);
    for(size_t j=r.begin(); j!=r.end();j++){
	
      double tmtp=mtp(j);
      if(tmtp>(-1)){
	tmtp=-1;
	mtp(j) =-1;
      }else{
	if(!std::isnan(tmtp)){
	  mu.col(j)=mu0.col(j)-2*tmtp*((mu1.col(j)-mu0.col(j)))+(tmtp*tmtp)*(mu.col(j)-2*mu1.col(j)+mu0.col(j));
	  //Rcpp::Rcout<<"SiRiSr_step:"<<std::endl;
	  SiRiSr.col(j)=(SiRiS.matrix()*(mu.col(j)).matrix()).array();
	}
      }
    }

    rassert(isnan(mu).sum()==0 && "mu has no NaN  values");
    rassert(isnan(SiRiSr).sum()==0 && "SiRiSr has no NaN  values");
  }
}

void rssr_norm::squarem_step_sparse(const blocked_range<size_t> & r) const{
    
  calculate_mtp(r);
  for(size_t j=r.begin(); j!=r.end();j++){
	
    double tmtp=mtp(j);
    if(tmtp>(-1)){
      tmtp=-1;
      mtp(j) =-1;
    }else{
      if(!std::isnan(tmtp)){
	mu.col(j)=mu0.col(j)-2*tmtp*((mu1.col(j)-mu0.col(j)))+(tmtp*tmtp)*(mu.col(j)-2*mu1.col(j)+mu0.col(j));
	//Rcpp::Rcout<<"SiRiSr_step:"<<std::endl;
	SiRiSr.col(j)=(sSiRiS*(mu.col(j)).matrix());
      }
    }
  }

  rassert(isnan(mu).sum()==0 && "mu has no NaN  values");
  rassert(isnan(SiRiSr).sum()==0 && "SiRiSr has no NaN  values");
}


      
void rssr_norm::squarem_backtrack(const blocked_range<size_t> & r)const {
  if(isSparse){
    squarem_backtrack_sparse(r);
  }else{
    for(size_t j=r.begin();j!=r.end();j++){
      double tmtp=mtp(j);
      double tlnZ=lnZ(j);
      double tlnZ0=lnZ0(j);
      if(tmtp<(-1) && (tlnZ < tlnZ0)){
	size_t num_bt=0;
	while((tlnZ<tlnZ0) && (num_bt < max_backtrack)){
	  tmtp=0.5*(tmtp-1);
	  mu.col(j) = mu0.col(j)-2*tmtp*(mu1.col(j)-mu0.col(j))+(tmtp*tmtp)*(mu.col(j)-2*mu1.col(j)+mu0.col(j));
	  SiRiSr.col(j) = (SiRiS*(mu.col(j).matrix())).array();
	  rss_vb_iter(blocked_range<size_t>(j,j),false);
	  rss_vb_iter(blocked_range<size_t>(j,j),true);
	  calc_lnZ(blocked_range<size_t>(j,j));
	  num_bt=num_bt+1;
	  final_btnum(j)++;
	}
	if(num_bt==max_backtrack){
	  mu.col(j)=mu0.col(j);
	  SiRiSr.col(j) = (SiRiS*(mu.col(j)).matrix()).array();
	  calc_lnZ(blocked_range<size_t>(j,j));
	  lnZ0(j)=lnZ00(j);
	}
      }
      mtp(j)=tmtp;
    }
  }
}

void rssr_norm::squarem_backtrack_sparse(const blocked_range<size_t> & r)const {
      
  for(size_t j=r.begin();j!=r.end();j++){
    double tmtp=mtp(j);
    double tlnZ=lnZ(j);
    double tlnZ0=lnZ0(j);
    if(tmtp<(-1) && (tlnZ < tlnZ0)){
      size_t num_bt=0;
      while((tlnZ<tlnZ0) && (num_bt < max_backtrack)){
	tmtp=0.5*(tmtp-1);
	mu.col(j) = mu0.col(j)-2*tmtp*(mu1.col(j)-mu0.col(j))+(tmtp*tmtp)*(mu.col(j)-2*mu1.col(j)+mu0.col(j));
	SiRiSr.col(j) = (sSiRiS*(mu.col(j).matrix()));
	rss_vb_iter(blocked_range<size_t>(j,j),false);
	rss_vb_iter(blocked_range<size_t>(j,j),true);
	calc_lnZ(blocked_range<size_t>(j,j));
	num_bt=num_bt+1;
	final_btnum(j)++;
      }
      if(num_bt==max_backtrack){
	mu.col(j)=mu0.col(j);
	SiRiSr.col(j) = (sSiRiS*(mu.col(j)).matrix());
	calc_lnZ(blocked_range<size_t>(j,j));
	lnZ0(j)=lnZ00(j);
      }
    }
    mtp(j)=tmtp;
  }
}
    
void rssr_norm::calc_err(const blocked_range<size_t> &r,double &max_err)const {
  Eigen::ArrayXd as(p);
  as.setOnes();
  max_err=0;
  if(use_lnztol){
	
    for(size_t j=r.begin(); j!=r.end();j++){
      double tmax_err=rel_err(lnZ(j),lnZ0(j));
      if(tmax_err>max_err){
	max_err=tmax_err;
      }
    }
  }else{
    for(size_t j=r.begin(); j!=r.end();j++){
      double tmax_err=find_maxerr(as,as,mu.col(j),mu0.col(j));
      if(tmax_err>max_err){
	max_err=tmax_err;
      }
    }
  }
}

void rssr_norm::compute_pve(const blocked_range<size_t> &r)const {

  for(size_t j=r.begin();j!=r.end();j++){
    pvevec(j)=(SiRiSr.col(j).matrix().transpose()*(mu.col(j)).matrix());
    pvevec(j)/=(double)n;
  }
} 



void rssr_norm::operator()(const blocked_range<size_t> &r)const {
    

  double max_err=1;
  size_t rstart=r.begin();
  size_t rsize= r.end()-rstart;
  calc_lnZ(r);
  int iter=0;

  
  ////Rcpp::Rcout<<"Initializing lnZ Segment "<<std::endl;    
  //  calc_lnZ(r);
  //<<"Initializing lnZ0 Segment "<<std::endl;
  lnZ0.segment(rstart,rsize)=lnZ.segment(rstart,rsize);
  while(max_err>tolerance){
    
    //<<"Setting lnZ00 Segment "<<std::endl;
    lnZ00.segment(rstart,rsize)=lnZ0.segment(rstart,rsize);
    //<<"Setting lnZ00 Segment "<<std::endl;
    lnZ0.segment(rstart,rsize)=lnZ.segment(rstart,rsize);
    //<<"Setting alpha0 block "<<std::endl;
    //<<"Setting mu0 block "<<std::endl;
    mu0.block(0,rstart,p,rsize)=mu.block(0,rstart,p,rsize);
    
    bool reverse = iter%2!=0;
    //<<"Performing first Iteration "<<std::endl;
    rss_vb_iter(r,reverse);
    //<<"Setting alpha1 block "<<std::endl;
    
    mu1.block(0,rstart,p,rsize)=mu.block(0,rstart,p,rsize);
    //<<"Performing second Iteration "<<std::endl;
    rss_vb_iter(r,reverse);
    
    
    rassert(isnan(mu).sum()==0 && "mu has no NaN  values");
    rassert(isnan(SiRiSr).sum()==0 && "SiRiSr has no NaN  values");
    
    
    
    rassert(isnan(mu0).sum()==0 && "mu0 has no NaN  values");
    rassert(isnan(mu1).sum()==0 && "mu1 has no NaN  values");
    
    
    
    rassert(isnan(lnZ).sum()==0 && "lnZ has no NaN  values");
    rassert(isnan(lnZ0).sum()==0 && "lnZ0 has no NaN  values");
    
    
    squarem_step(r);
    //<<"Third Iteration "<<std::endl;
    rss_vb_iter(r,reverse);
    
    calc_lnZ(r);
    //<<"Squarem Backtrack "<<std::endl;
    squarem_backtrack(r);
    //<<"Error Calculation "<<std::endl;
    calc_err(r,max_err);
    if(iter>itermax){
      for(size_t  j=r.begin();j!=r.end();j++){
	final_max_err(j)=max_err;
	iternum(j)=iter;
      }
      compute_pve(r);
      break;
    }
    iter++;
  }
  
  for(size_t  j=r.begin();j!=r.end();j++){
    final_max_err(j)=max_err;
    iternum(j)=iter;
    alpha_mean(j)=1;
    mu_mean(j)=mu.col(j).mean();
  }
  rassert(isnan(alpha_mean).sum()==0 && "alpha_mean is entirely finite");
  compute_pve(r);
}
