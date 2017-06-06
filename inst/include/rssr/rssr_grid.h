#ifndef RSSR_GRID_H
#define RSSR_GRID_H

#include <RcppEigen.h>
#include "rssr_types.h"



double initialize_value(const double val, const double param);

double calculate_mtp(const c_arrayxd_internal alpha,const c_arrayxd_internal alpha1,const c_arrayxd_internal alpha0, const c_arrayxd_internal mu,const c_arrayxd_internal mu1, const c_arrayxd_internal mu0);

double calculate_mtp(const c_arrayxd_internal alpha_r,const c_arrayxd_internal alpha_v,const c_arrayxd_internal mu_r, const c_arrayxd_internal mu_v);

Eigen::ArrayXd initialize_array(const c_arrayxd_internal init, const double param);

void calc_max_err(const double lnZ,const double lnZ0,const c_arrayxd_internal alpha, const c_arrayxd_internal alpha0,const c_arrayxd_internal mu,const c_arrayxd_internal mu0,  double &max_err,const bool lnztol);

Eigen::ArrayXd initialize_s(const c_arrayxd_internal sesquare, const double sigma_beta_square);



Eigen::ArrayXd initialize_ssrat(const c_arrayxd_internal s, const double sigma_beta_square);

Rcpp::DataFrame grid_rss_varbvsr(
    const c_sparseMatrix_internal  SiRiS,
    const c_arrayxd_internal sigma_beta,
    const c_arrayxd_internal logodds,
    const c_arrayxd_internal  betahat,
    const c_arrayxd_internal se,
    const c_arrayxd_internal talpha0,
    const c_arrayxd_internal tmu0,
    const c_arrayxd_internal tSiRiSr0,
    double tolerance,
    int itermax,
    bool verbose,
    bool lnz_tol);

Rcpp::DataFrame grid_search_rss_varbvsr_sp(
    const sparseMatrix_external SiRiS,
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
    Rcpp::LogicalVector lnz_tol);


Rcpp::DataFrame grid_rss_varbvsr(
    const c_Matrix_internal  SiRiS,
    const c_arrayxd_internal sigma_beta,
    const c_arrayxd_internal logodds,
    const c_arrayxd_internal  betahat,
    const c_arrayxd_internal se,
    const c_arrayxd_internal talpha0,
    const c_arrayxd_internal tmu0,
    const c_arrayxd_internal tSiRiSr0,
    double tolerance,
    int itermax,
    bool verbose,
    bool lnz_tol);

Rcpp::DataFrame grid_search_rss_varbvsr(
    const Matrix_external SiRiS,
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
    Rcpp::LogicalVector lnz_tol);

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
                                              Rcpp::LogicalVector lnz_tol);


Rcpp::DataFrame grid_rss_varbvsr_array(
    const c_Matrix_internal SiRiS,
    const c_arrayxd_internal sigma_beta,
    const c_arrayxd_internal logodds,
    const c_arrayxd_internal  betahat,
    const c_arrayxd_internal  se,
    const c_arrayxd_internal talpha0,
    const c_arrayxd_internal tmu0,
    const c_arrayxd_internal tSiRiSr0,
    double tolerance,
    int itermax,
    bool isVerbose,
    bool islnz_tol);


Rcpp::DataFrame grid_rss_varbvsr_serial(
    const c_Matrix_internal SiRiS,
    const c_arrayxd_internal sigma_beta,
    const c_arrayxd_internal logodds,
    const c_arrayxd_internal  betahat,
    const c_arrayxd_internal  se,
    const c_arrayxd_internal talpha0,
    const c_arrayxd_internal tmu0,
    const c_arrayxd_internal tSiRiSr0,
    double tolerance,
    int itermax,
    bool isVerbose,
    bool islnz_tol);


#endif
