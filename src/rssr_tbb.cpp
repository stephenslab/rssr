#include <RcppEigen.h>
#include "rssr.h"

//[[Rcpp::depends(RcppParallel)]]
using namespace Rcpp;
using namespace RcppParallel;
using namespace tbb;


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


class EigenVectord{
public:
    const size_t p;
  Rcpp::NumericVector rvec;
  RcppParallel::RVector<double> pvec;
  mdarray evec;

  EigenVectord(const size_t p_):p(p_),
				rvec(p),
				pvec(rvec),
				evec(rvec.begin(),p)
  {};
};

class EigenVectori{
public:
    const size_t p;
  Rcpp::IntegerVector rvec;
  RcppParallel::RVector<int> pvec;
  miarray evec;

  EigenVectori(const size_t p_):p(p_),
				rvec(p),
				pvec(rvec),
				evec(rvec.begin(),p)
  {};
};

class ERvectord{
public:
    const size_t p;
  RcppParallel::RVector<double> pvec;
  mdarray evec;

  ERvectord(NumericVector &rvec,const size_t p_):p(p_),
							  pvec(rvec),
							  evec(rvec.begin(),p)
  {
    rassert(evec(0)==rvec[0] &&" evec and rvec are the same!");
  };
};

class ERvectori{
public:
    const size_t p;
  RcppParallel::RVector<int> pvec;
  miarray evec;

  ERvectori(IntegerVector &rvec,const size_t p_):p(p_),
							  pvec(rvec),
							  evec(rvec.begin(),p)
  {
    rassert(evec(0)==rvec[0] &&" evec and rvec are the same!");
  };
};

class cERvectord{
public:
  const size_t p;
  RcppParallel::RVector<double> pvec;
  const c_mdarray evec;

  cERvectord(const NumericVector &rvec,const size_t p_):p(p_),
							      pvec(rvec),
							      evec(rvec.begin(),p)
  {
    rassert(evec(0)==rvec[0] &&" evec and rvec are the same!");
  };
};

class EigenMatrixd{
public:
    const size_t rows;
  const size_t cols;
  Rcpp::NumericMatrix rmat;
  RcppParallel::RMatrix<double> pmat;
  m2darray emat;

  EigenMatrixd(const size_t rows_,const size_t cols_):rows(rows_),
						      cols(cols_),
						      rmat(rows,cols),
						      pmat(rmat),
						      emat(rmat.begin(),rows,cols)
  {    rassert(emat(0,0)==rmat(0,0) &&" emat and rmat are the same!");};
};
    
class ERmatrixd{
public:
    const size_t rows;
  const size_t cols;
  RcppParallel::RMatrix<double> pmat;
  m2darray emat;

  ERmatrixd(NumericMatrix &rmat,const size_t rows_,const size_t cols_):rows(rows_),
									       cols(cols_),
					     
									       pmat(rmat),
									       emat(rmat.begin(),rows,cols)
  {};
};

class cERmatrixd{
public:
    const size_t rows;
  const size_t cols;
  RcppParallel::RMatrix<double> pmat;
  const c_mmat emat;

  cERmatrixd(const NumericMatrix &rmat,const size_t rows_,const size_t cols_):rows(rows_),
									     cols(cols_),
									     pmat(rmat),
									     emat(rmat.begin(),rows,cols)
  {};
};

class cERarrayd{
public:
  const size_t rows;
  const size_t cols;
  RcppParallel::RMatrix<double> pmat;
  const c_m2darray earray;

  cERarrayd(const NumericMatrix &rmat,const size_t rows_,const size_t cols_):rows(rows_),
									     cols(cols_),
									     pmat(rmat),
									     earray(rmat.begin(),rows,cols)
  {};
};

class ERarrayd{
public:
    const size_t rows;
  const size_t cols;
  RcppParallel::RMatrix<double> pmat;
  const m2darray earray;

  ERarrayd(NumericMatrix &rmat,const size_t rows_,const size_t cols_):rows(rows_),
									     cols(cols_),
									     pmat(rmat),
									     earray(rmat.begin(),rows,cols)
  {};
};




template<typename T, typename D> void zero_out(T &dyn_array, const D fillval){
#ifndef EIGEN_NO_DEBUG
  std::fill(dyn_array.local().begin(),dyn_array.local().end(),fillval);
#endif
}
  

using namespace Eigen;




template<typename T> class Data_wrap{
public:
  c_mdarray betahat;
  c_mdarray se;
  mdarray sesquare;
  mdarray q;
  const int n;
  const T SiRiS;
  
    
  Data_wrap(c_mdarray betahat_,
	    c_mdarray se_,
	    mdarray sesquare_,
	    mdarray q_,
	    const T SiRiS_,
	    const int n_): betahat(betahat_),
													    se(se_),
													    sesquare(sesquare_),
													    q(q_),
													    SiRiS(SiRiS_),n(n_){};
};

class Param_wrap{
public:
  c_mdarray logodds;
  c_mdarray sigma_beta;
  mdarray sigma_beta_square;
  c_m2darray s;
  c_m2darray ssrat;
  std::vector<int> *forward_range;
  std::vector<int> *reverse_range;
  
  Param_wrap(c_mdarray logodds_,
	     c_mdarray sigma_beta_,
	     mdarray sigma_beta_square_,
	     c_m2darray s_,
	     c_m2darray ssrat_,
	     std::vector<int> *forward_range_,
	     std::vector<int> *reverse_range_): logodds(logodds_),
											     sigma_beta(sigma_beta_),
											     sigma_beta_square(sigma_beta_square_),
											     s(s_),
											     ssrat(ssrat_),forward_range(forward_range_),
											     reverse_range(reverse_range_){};
};
  
  
  
class Fit_res{
 public:
  const size_t gridsize;
  const size_t p;
  mdarray lnZ;
  mdarray lnZ0;
  mdarray lnZ00;
  mdarray alpha_mean;
  mdarray mu_mean;
  mdarray pvevec;
  miarray iternum;
  mdarray final_max_err;
  miarray final_btnum;
  mdarray mtp;
  
  Fit_res(mdarray lnZ_,
	  mdarray lnZ0_,
	  mdarray lnZ00_,
	  mdarray alpha_mean_,
	  mdarray mu_mean_,
	  mdarray pvevec_,
	  miarray  iternum_,
	  mdarray final_max_err_,
	  miarray final_btnum_,
	  mdarray mtp_,const size_t p_,const size_t gridsize_):  lnZ(lnZ_),
								 lnZ0(lnZ0_),
								 lnZ00(lnZ00_),
								 alpha_mean(alpha_mean_),
								 mu_mean(mu_mean_),
								 pvevec(pvevec_),
								 iternum(iternum_),
								 final_max_err(final_max_err_),
								 final_btnum(final_btnum_),
								 mtp(mtp_),
								 p(p_),
								 gridsize(gridsize_){};

};


class Temp_fit{
public:

  const size_t gridsize;
  const size_t p;
  
  m2darray alpha;
  m2darray alpha0;
  m2darray alpha1;

  m2darray mu;
  m2darray mu0;
  m2darray mu1;
  mdarray rvec;
  m2darray SiRiSr;

  Temp_fit(c_mdarray oalpha,
	   c_mdarray omu,
	   c_mdarray oSiRiSr,
	   m2darray alpha_,
	   m2darray alpha0_,
	   m2darray alpha1_,
	   m2darray mu_,
	   m2darray mu0_,
	   m2darray mu1_,
	   m2darray SiRiSr_,
	   mdarray rvec_,
	   const size_t p_,size_t gridsize_):gridsize(gridsize_),
					     p(p_),
					     alpha(alpha_),
					     mu(mu_),
					     SiRiSr(SiRiSr_),
					     alpha0(alpha0_),
					     alpha1(alpha1_),
					     mu0(mu0_),
					     mu1(mu1_),
					     rvec(rvec_)

  {
    for(size_t i=0;i<gridsize;i++){
      alpha.col(i)=oalpha;
      mu.col(i)=omu;
      SiRiSr.col(i)=oSiRiSr;
    }
  };

  template<typename T> void calc_lnZ(const blocked_range<size_t> & r, const Data_wrap<T> &data_wrap,
				     const Param_wrap &param_wrap,mdarray &lnZ){
    for(size_t j=r.begin(); j!=r.end();j++){

      
      double tlogodds=param_wrap.logodds(j);


      lnZ(j)=calculate_lnZ(data_wrap.q,
			   alpha.col(j)*mu.col(j),
			   SiRiSr.col(j),
			   tlogodds,
			   data_wrap.sesquare,
			   alpha.col(j),
			   mu.col(j),
			   param_wrap.s.col(j),
			   param_wrap.sigma_beta(j));
      
      if(std::isnan(lnZ(j))){
	std::cerr<<j<<std::endl;
	std::cerr<<"lnZ(j):"<<lnZ(j)<<std::endl;
	rassert(!isnan(lnZ(j))&& "lnZ(j) is finite ");
      }
    }
    };
  void calculate_mtp(const blocked_range <size_t > &r,mdarray &mtp){
    
    for(size_t i=r.begin();i!=r.end();i++){
      mtp(i)=
	-std::sqrt(((alpha1.col(i)-alpha0.col(i)).square()).sum()+((mu1.col(i)-mu0.col(i)).square()).sum())/std::sqrt(((alpha.col(i)-2*alpha1+alpha0).square()).sum()+((mu.col(i)-2*mu1.col(i)+mu0.col(i)).square()).sum()+double_lim::epsilon());
    }
  };
  template <typename T> void rss_vb_iter(const blocked_range<size_t> r,bool reverse, const Param_wrap &param_wrap, const Data_wrap<T> &data_wrap){
    size_t rstart=r.begin();
    size_t rsize= r.end()-rstart;

    mdarray q(data_wrap.q);
    c_m2darray s(param_wrap.s);
    c_m2darray ssrat(param_wrap.ssrat);
    T SiRiS(data_wrap.SiRiS);
    mdarray sesquare(data_wrap.sesquare);
    c_mdarray logodds(param_wrap.logodds);
    std::vector<int> *range;
    if(reverse){
      range=param_wrap.reverse_range;
    }else{
      range =param_wrap.forward_range;
    }
    for(const int &i: *range){
      for(size_t j=r.begin();j!=r.end();j++){
	rvec(j)=alpha(i,j)*mu(i,j);
	mu(i,j)=s(i,j)*(q(i)+rvec(j)/sesquare(i)-SiRiSr(i,j));
	alpha(i,j)=sigmoid(logodds(j)+0.5*(ssrat(i,j)+(mu(i,j)*mu(i,j))/s(i,j)));
	SiRiSr.matrix().col(j)+=SiRiS.col(i)*(alpha(i,j)*mu(i,j)-rvec(j));
      }
    }
  };
  template<typename T> void squarem_step(const blocked_range<size_t> & r,mdarray &mtp, const T &SiRiS){
    calculate_mtp(r,mtp);
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
	  SiRiSr.col(j)=(SiRiS*(alpha.col(j)*mu.col(j)).matrix()).array();
	}
      }
    }
    rassert(isnan(alpha).sum()==0 && "alpha has no NaN  values");
    rassert(isnan(mu).sum()==0 && "mu has no NaN  values");
    rassert(isnan(SiRiSr).sum()==0 && "SiRiSr has no NaN  values");
  }

  template<typename T> void squarem_backtrack(const blocked_range<size_t> & r, const Data_wrap<T> &data_wrap,Fit_res &fit_res, const Param_wrap &param_wrap, const size_t max_backtrack=10){
    for(size_t j=r.begin();j!=r.end();j++){
      double tmtp=fit_res.mtp(j);
      double tlnZ=fit_res.lnZ(j);
      double tlnZ0=fit_res.lnZ0(j);
      const T SiRiS(data_wrap.SiRiS);
      if(tmtp<(-1) && (tlnZ < tlnZ0)){
	size_t num_bt=0;
	while((tlnZ<tlnZ0) && (num_bt < max_backtrack)){
	  tmtp=0.5*(tmtp-1);
	  alpha.col(j) = alpha0.col(j)-2*tmtp*(alpha1.col(j)-alpha0.col(j))+(tmtp*tmtp)*(alpha.col(j)-2*alpha1.col(j)+alpha0.col(j));
	  mu.col(j) = mu0.col(j)-2*tmtp*(mu1.col(j)-mu0.col(j))+(tmtp*tmtp)*(mu.col(j)-2*mu1.col(j)+mu0.col(j));
	  SiRiSr.col(j) = (SiRiS*(alpha.col(j)*mu.col(j)).matrix()).array();
	  rss_vb_iter(blocked_range<size_t>(j,j),false,param_wrap,data_wrap);
	  rss_vb_iter(blocked_range<size_t>(j,j),true,param_wrap,data_wrap);
	  
	  calc_lnZ(blocked_range<size_t>(j,j),data_wrap,param_wrap,fit_res.lnZ);
	  num_bt=num_bt+1;
	  fit_res.final_btnum(j)+=1;
	}
	if(num_bt==max_backtrack){
	  alpha.col(j)=alpha0.col(j);
	  mu.col(j)=mu0.col(j);
	  SiRiSr.col(j) = (SiRiS*(alpha.col(j)*mu.col(j)).matrix()).array();
	  calc_lnZ(blocked_range<size_t>(j,j),data_wrap,param_wrap,fit_res.lnZ);
	  fit_res.lnZ0(j)=fit_res.lnZ00(j);
	}
      }
      fit_res.mtp(j)=tmtp;
    }
  }
  void compute_pve(const blocked_range<size_t> &r,mdarray &pvevec,const int n=1){
    for(size_t j=r.begin();j!=r.end();j++){
      pvevec(j)=(SiRiSr.col(j).matrix().transpose()*(alpha.col(j)*mu.col(j)).matrix());
      pvevec(j)/=(double)n;
    }
  };    
};


template <typename T>
struct is_sparse{
  static const bool value = false;
};
template <>
struct is_sparse<Eigen::Map< Eigen::SparseMatrix<double> > >{
  static const bool value = true;
};

	   
template <typename T> class RSSR{
public:
  mutable Temp_fit temp_fit;
  mutable Fit_res fit_res;
  mutable Param_wrap param_wrap;
  mutable Data_wrap<T> data_wrap;
  const bool doAccelerate;
  const bool isSparse;
  const bool use_lnztol;
  const double tolerance;
  const size_t p;
  const size_t gridsize;
  const size_t itermax;
  RSSR(Temp_fit &temp_fit_,
       Fit_res &fit_res_,
       Param_wrap &param_wrap_,
       Data_wrap<T> &data_wrap_,
       const bool doAccelerate_,
       const bool use_lnztol_=true,
       const double tolerance_=0.001,const size_t itermax_=100): temp_fit(temp_fit_),
								 fit_res(fit_res_),
								 param_wrap(param_wrap_),
								 data_wrap(data_wrap_),
								 doAccelerate(doAccelerate_),
								 isSparse(is_sparse<T>::value),
								 use_lnztol(use_lnztol_),
								 tolerance(tolerance_),
								 gridsize(temp_fit.alpha.cols()),
								 itermax(itermax_),
								  p(temp_fit.alpha.rows()){};

  void calc_err(const blocked_range<size_t> &r,double &max_err)const {
  
    max_err=0;
    if(use_lnztol){
      for(size_t j=r.begin(); j!=r.end();j++){
	double tmax_err=rel_err(fit_res.lnZ(j),fit_res.lnZ0(j));
	if(tmax_err>max_err){
	  max_err=tmax_err;
	}
      }
    }else{
      for(size_t j=r.begin(); j!=r.end();j++){
	double tmax_err=find_maxerr(temp_fit.alpha.col(j),temp_fit.alpha0.col(j),temp_fit.alpha.col(j)*temp_fit.mu.col(j),temp_fit.alpha0.col(j)*temp_fit.mu0.col(j));
	if(tmax_err>max_err){
	  max_err=tmax_err;
	}
      }
    }
  };
  void operator()(    const blocked_range<size_t> &r)const{    
    double max_err=1;
    size_t rstart=r.begin();
    size_t rsize= r.end()-rstart;
    //  blocked_range<size_t> r(0,rsize);
    //  tls_init(r);
    temp_fit.calc_lnZ(r,data_wrap,param_wrap,fit_res.lnZ);
    int iter=0;
  
    // Rcpp::Rcout<<"Initializing lnZ Segment "<<std::endl;    
    //  calc_lnZ(r);
    //Rcpp::Rcout<<"Initializing lnZ0 Segment "<<std::endl;
    fit_res.lnZ0.segment(rstart,rsize)=fit_res.lnZ.segment(rstart,rsize);
    while(max_err>tolerance){
    
      //Rcpp::Rcout<<"Setting lnZ00 Segment "<<std::endl;
      fit_res.lnZ00.segment(rstart,rsize)=fit_res.lnZ0.segment(rstart,rsize);
      //Rcpp::Rcout<<"Setting lnZ00 Segment "<<std::endl;
      fit_res.lnZ0.segment(rstart,rsize)=fit_res.lnZ.segment(rstart,rsize);
      //Rcpp::Rcout<<"Setting alpha0 block "<<std::endl;
      temp_fit.alpha0.block(0,rstart,p,rsize)=temp_fit.alpha.block(0,rstart,p,rsize);
      //Rcpp::Rcout<<"Setting mu0 block "<<std::endl;
      temp_fit.mu0.block(0,rstart,p,rsize)=temp_fit.mu.block(0,rstart,p,rsize);
    
      bool reverse = iter%2!=0;
      //Rcpp::Rcout<<"Performing first Iteration "<<std::endl;
      temp_fit.rss_vb_iter(r,reverse,param_wrap,data_wrap);
      //Rcpp::Rcout<<"Setting alpha1 block "<<std::endl;
      temp_fit.alpha1.block(0,rstart,p,rsize)=temp_fit.alpha.block(0,rstart,p,rsize);
      temp_fit.mu1.block(0,rstart,p,rsize)=temp_fit.mu.block(0,rstart,p,rsize);
      //Rcpp::Rcout<<"Performing second Iteration "<<std::endl;
      temp_fit.rss_vb_iter(r,reverse,param_wrap,data_wrap);
    

      //Rcpp::Rcout<<"Squarem Step "<<std::endl;
      temp_fit.squarem_step(r,fit_res.mtp,data_wrap.SiRiS);
      //Rcpp::Rcout<<"Third Iteration "<<std::endl;
      temp_fit.rss_vb_iter(r,reverse,param_wrap,data_wrap);
      temp_fit.calc_lnZ(r,data_wrap,param_wrap,fit_res.lnZ);
      //Rcpp::Rcout<<"Squarem Backtrack "<<std::endl;
      temp_fit.squarem_backtrack(r,data_wrap,fit_res,param_wrap,10);
      //Rcpp::Rcout<<"Error Calculation "<<std::endl;
      calc_err(r,max_err);
      if(iter>itermax){
	for(size_t  j=r.begin();j!=r.end();j++){
	  fit_res.final_max_err(j)=max_err;
	  fit_res.iternum(j)=iter;
	}
	temp_fit.compute_pve(r,fit_res.pvevec,data_wrap.n);
	//      results->write_res(qr);
	break;
      
      }
      iter++;
    }
    //Rcpp::Rcout<<"Error Calculation "<<std::endl;
    
    for(size_t  j=r.begin();j!=r.end();j++){
      fit_res.final_max_err(j)=max_err;
      fit_res.iternum(j)=iter;
      fit_res.alpha_mean(j)=temp_fit.alpha.col(j).mean();
      fit_res.mu_mean(j)=temp_fit.mu.col(j).mean(); 
    }


    temp_fit.compute_pve(r,fit_res.pvevec,data_wrap.n);
    //  results->write_res(qr);
  };
};

  
    
  








template <typename T>Rcpp::DataFrame grid_search_rss_varbvsr(
							     const  T SiRiS,
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

  

  cERvectord tsigma_beta(sigma_beta,tot_size);
  cERvectord tlogodds(logodds,tot_size);
  cERvectord tbetahat(betahat,p);
  cERvectord tse(se,p);
  cERvectord talpha0(alpha0,p);
  cERvectord tmu0(mu0,p);
  cERvectord tSiRiSr0(SiRiSr0,p);

  //std::cout<<"allocating tsigma_beta_square"<<std::endl;
  NumericVector tsigma_beta_square(tot_size);
    //std::cout<<"mapping sigma_beta_square"<<std::endl;
  ERvectord sigma_beta_square(tsigma_beta_square,tot_size);
  //std::cout<<"setting sigma_beta_square"<<std::endl;
  //  std::cout<<tsigma_beta.evec<<std::endl;
  //  std::cout<<"sigma_beta!"<<std::endl;
  //  std::cout<<sigma_beta_square.evec<<std::endl;
  
  sigma_beta_square.evec=tsigma_beta.evec*tsigma_beta.evec;
  //std::cout<<"allocating tsesquare"<<std::endl;
  NumericVector tsesquare(p);
  ERvectord sesquare(tsesquare,p);
  if(sesquare.evec(0)!=tsesquare(0)){
    Rcpp::stop("Something is weird with the ERvectord class");
  }
  sesquare.evec=tse.evec*tse.evec;

  NumericVector tq(p);
  ERvectord q(tq,p);
  q.evec=tbetahat.evec/sesquare.evec;
  
  //std::cout<<"Initializing s"<<std::endl;
  NumericMatrix s=initialize_smat(tsesquare,tsigma_beta_square);
  //std::cout<<"Initializing ssrat"<<std::endl;
  NumericMatrix ssrat=initialize_ssratmat(s,tsigma_beta_square);

  cERarrayd ts(s,p,tot_size);
  cERarrayd tssrat(ssrat,p,tot_size);


  NumericMatrix talphas(p,tot_size);
  NumericMatrix talpha0s(p,tot_size);
  NumericMatrix talpha1s(p,tot_size);
  
  ERarrayd alphas(talphas,p,tot_size);
  ERarrayd alpha0s(talpha0s,p,tot_size);
  ERarrayd alpha1s(talpha1s,p,tot_size);


  NumericMatrix tmus(p,tot_size);
  NumericMatrix tmu0s(p,tot_size);
  NumericMatrix tmu1s(p,tot_size);
  
  ERarrayd mus(tmus,p,tot_size);
  ERarrayd mu0s(tmu0s,p,tot_size);
  ERarrayd mu1s(tmu1s,p,tot_size);


  NumericMatrix tSiRiSrs(p,tot_size);
  ERarrayd SiRiSrs(tSiRiSrs,p,tot_size);
  NumericVector trvecs(tot_size);
  ERvectord rvecs(trvecs,tot_size);
  


  Temp_fit temp_fit(talpha0.evec,
		    tmu0.evec,
		    tSiRiSr0.evec,
		    alphas.earray,
		    alpha0s.earray,
		    alpha1s.earray,
		    mus.earray,
		    mu0s.earray,
		    mu1s.earray,
		    SiRiSrs.earray,
		    rvecs.evec,p,tot_size);
		    



  

  IntegerVector rfinal_btnum(tot_size);
  IntegerVector riternum(tot_size);
  NumericVector ralpha_mean(tot_size);
  NumericVector rmu_mean(tot_size);
  NumericVector rpvevec(tot_size);
  NumericVector rlnZ(tot_size);
  NumericVector rlnZ0(tot_size);
  NumericVector rlnZ00(tot_size);
  NumericVector rmtp(tot_size);
  NumericVector rfinal_max_err(tot_size);
  
  ERvectori trfinal_btnum(rfinal_btnum,tot_size);
  ERvectori triternum(riternum,tot_size);
  ERvectord tralpha_mean(ralpha_mean,tot_size);
  ERvectord trmu_mean(rmu_mean,tot_size);
  ERvectord trpvevec(rpvevec,tot_size);
  ERvectord trlnZ(rlnZ,tot_size);
  ERvectord trlnZ0(rlnZ0,tot_size);
  ERvectord trlnZ00(rlnZ00,tot_size);
  ERvectord trmtp(rmtp,tot_size);
  ERvectord trfinal_max_err(rfinal_max_err,tot_size);
  



  std::vector<int> reverse_range(p);
  std::vector<int> forward_range(p);
  for(size_t i=0;i<p;i++){
    forward_range[i]=i;
  }
  std::reverse_copy(forward_range.begin(),forward_range.end(),reverse_range.begin());  

  
  Data_wrap<T> data_wrap(tbetahat.evec,
			 tse.evec,
			 sesquare.evec,
			 q.evec,
			 SiRiS,
			 n);
  Param_wrap param_wrap(tlogodds.evec,
			tsigma_beta.evec,
			sigma_beta_square.evec,
			ts.earray,
			tssrat.earray,
			&forward_range,
			&reverse_range);
  Fit_res fit_res(trlnZ.evec,
		  trlnZ0.evec,
		  trlnZ00.evec,
		  tralpha_mean.evec,
		  trmu_mean.evec,
		  trpvevec.evec,
		  triternum.evec,
		  trfinal_max_err.evec,
		  trfinal_btnum.evec,
		  trmtp.evec,
		  p,
		  tot_size);
  

  


  //std::cout  <<"Initializing RSSR"<<std::endl;
  RSSR<T> rssr_obj(temp_fit,
		   fit_res,
		   param_wrap,
		   data_wrap,
		   true,islnz_tol,tolerance);
  //Rcpp::Rcout  <<"Starting parallel Algorithm!"<<std::endl;
  //rssr_obj(0,tot_size);
  //  Rcpp::Rcout<<"Starting Serial Algorithm!"<<std::endl;
  //rssr_obj(blocked_range<size_t>(0,tot_size));
  parallel_for(blocked_range<size_t>(0,tot_size,grainsize),rssr_obj);
  //    Rcpp::Rcout<<"Finished!"<<std::endl;
  //  std::cout<<"Checking then Returning!"<<std::endl;

    //std::cout<<"Returning!"<<std::endl;
  
    
  return(Rcpp::DataFrame::create(_["logodds"]=logodds,
				 _["sigb"]=sigma_beta,
				 _["rel_err"]=rfinal_max_err,
				 _["iterations"]=riternum,
				 _["alpha_mean"]=ralpha_mean,
				 _["mu_mean"]=rmu_mean,
				 _["pve"]=rpvevec,
				 _["lnZ"]=rlnZ));
}




//[[Rcpp::export(name="grid_search_rss_varbvsr_dense")]]
Rcpp::DataFrame grid_search_rss_varbvsr_d(
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
  
  size_t p=betahat.size();
  size_t tot_size=logodds.size();
  cERmatrixd tSiRiS(SiRiS,p,p);

  return(grid_search_rss_varbvsr<c_mmat>(tSiRiS.emat,
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


//[[Rcpp::export(name="grid_search_rss_varbvsr_sparse")]]
Rcpp::DataFrame grid_search_rss_varbvsr_s(
    const  sparseMatrix_external &SiRiS,
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
  
  size_t p=betahat.size();
  size_t tot_size=logodds.size();

  return(grid_search_rss_varbvsr<sparseMatrix_external>(SiRiS,
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
