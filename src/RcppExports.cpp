// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/rssr.h"
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// logsigmoid_exp
Eigen::ArrayXd logsigmoid_exp(const arrayxd_external x);
RcppExport SEXP rssr_logsigmoid_exp(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arrayxd_external >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(logsigmoid_exp(x));
    return rcpp_result_gen;
END_RCPP
}
// sigmoid_exp
Eigen::ArrayXd sigmoid_exp(const arrayxd_external x);
RcppExport SEXP rssr_sigmoid_exp(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arrayxd_external >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(sigmoid_exp(x));
    return rcpp_result_gen;
END_RCPP
}
// exp_betavar
Eigen::ArrayXd exp_betavar(const arrayxd_external p, const arrayxd_external mu, const arrayxd_external s);
RcppExport SEXP rssr_exp_betavar(SEXP pSEXP, SEXP muSEXP, SEXP sSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arrayxd_external >::type p(pSEXP);
    Rcpp::traits::input_parameter< const arrayxd_external >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const arrayxd_external >::type s(sSEXP);
    rcpp_result_gen = Rcpp::wrap(exp_betavar(p, mu, s));
    return rcpp_result_gen;
END_RCPP
}
// exp_intklbeta_rssbvsr
double exp_intklbeta_rssbvsr(const arrayxd_external alpha, const arrayxd_external mu, const arrayxd_external sigma_square, double sigma_beta_square);
RcppExport SEXP rssr_exp_intklbeta_rssbvsr(SEXP alphaSEXP, SEXP muSEXP, SEXP sigma_squareSEXP, SEXP sigma_beta_squareSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arrayxd_external >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const arrayxd_external >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const arrayxd_external >::type sigma_square(sigma_squareSEXP);
    Rcpp::traits::input_parameter< double >::type sigma_beta_square(sigma_beta_squareSEXP);
    rcpp_result_gen = Rcpp::wrap(exp_intklbeta_rssbvsr(alpha, mu, sigma_square, sigma_beta_square));
    return rcpp_result_gen;
END_RCPP
}
// exp_intgamma
double exp_intgamma(double logodds, const arrayxd_external alpha);
RcppExport SEXP rssr_exp_intgamma(SEXP logoddsSEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type logodds(logoddsSEXP);
    Rcpp::traits::input_parameter< const arrayxd_external >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(exp_intgamma(logodds, alpha));
    return rcpp_result_gen;
END_RCPP
}
// rel_err
double rel_err(double p0, double p1);
RcppExport SEXP rssr_rel_err(SEXP p0SEXP, SEXP p1SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type p0(p0SEXP);
    Rcpp::traits::input_parameter< double >::type p1(p1SEXP);
    rcpp_result_gen = Rcpp::wrap(rel_err(p0, p1));
    return rcpp_result_gen;
END_RCPP
}
// exp_find_maxerr
double exp_find_maxerr(const arrayxd_external alpha, const arrayxd_external alpha0, const arrayxd_external r, const arrayxd_external r0);
RcppExport SEXP rssr_exp_find_maxerr(SEXP alphaSEXP, SEXP alpha0SEXP, SEXP rSEXP, SEXP r0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arrayxd_external >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const arrayxd_external >::type alpha0(alpha0SEXP);
    Rcpp::traits::input_parameter< const arrayxd_external >::type r(rSEXP);
    Rcpp::traits::input_parameter< const arrayxd_external >::type r0(r0SEXP);
    rcpp_result_gen = Rcpp::wrap(exp_find_maxerr(alpha, alpha0, r, r0));
    return rcpp_result_gen;
END_RCPP
}
// exp_update_logodds
double exp_update_logodds(const arrayxd_external alpha);
RcppExport SEXP rssr_exp_update_logodds(SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arrayxd_external >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(exp_update_logodds(alpha));
    return rcpp_result_gen;
END_RCPP
}
// exp_calculate_lnZ
double exp_calculate_lnZ(const vectorxd_external q, const vectorxd_external r, const vectorxd_external SiRiSr, double logodds, const vectorxd_external sesquare, const vectorxd_external alpha, const vectorxd_external mu, const vectorxd_external s, double sigb);
RcppExport SEXP rssr_exp_calculate_lnZ(SEXP qSEXP, SEXP rSEXP, SEXP SiRiSrSEXP, SEXP logoddsSEXP, SEXP sesquareSEXP, SEXP alphaSEXP, SEXP muSEXP, SEXP sSEXP, SEXP sigbSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const vectorxd_external >::type q(qSEXP);
    Rcpp::traits::input_parameter< const vectorxd_external >::type r(rSEXP);
    Rcpp::traits::input_parameter< const vectorxd_external >::type SiRiSr(SiRiSrSEXP);
    Rcpp::traits::input_parameter< double >::type logodds(logoddsSEXP);
    Rcpp::traits::input_parameter< const vectorxd_external >::type sesquare(sesquareSEXP);
    Rcpp::traits::input_parameter< const vectorxd_external >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const vectorxd_external >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const vectorxd_external >::type s(sSEXP);
    Rcpp::traits::input_parameter< double >::type sigb(sigbSEXP);
    rcpp_result_gen = Rcpp::wrap(exp_calculate_lnZ(q, r, SiRiSr, logodds, sesquare, alpha, mu, s, sigb));
    return rcpp_result_gen;
END_RCPP
}
// rss_varbvsr_naive
Rcpp::List rss_varbvsr_naive(const Matrix_external SiRiS, const double sigma_beta, const double logodds, const arrayxd_external betahat, const arrayxd_external se, const arrayxd_external talpha0, const arrayxd_external tmu0, const arrayxd_external tSiRiSr0, double tolerance, int itermax, Rcpp::LogicalVector verbose, Rcpp::LogicalVector lnz_tol);
RcppExport SEXP rssr_rss_varbvsr_naive(SEXP SiRiSSEXP, SEXP sigma_betaSEXP, SEXP logoddsSEXP, SEXP betahatSEXP, SEXP seSEXP, SEXP talpha0SEXP, SEXP tmu0SEXP, SEXP tSiRiSr0SEXP, SEXP toleranceSEXP, SEXP itermaxSEXP, SEXP verboseSEXP, SEXP lnz_tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Matrix_external >::type SiRiS(SiRiSSEXP);
    Rcpp::traits::input_parameter< const double >::type sigma_beta(sigma_betaSEXP);
    Rcpp::traits::input_parameter< const double >::type logodds(logoddsSEXP);
    Rcpp::traits::input_parameter< const arrayxd_external >::type betahat(betahatSEXP);
    Rcpp::traits::input_parameter< const arrayxd_external >::type se(seSEXP);
    Rcpp::traits::input_parameter< const arrayxd_external >::type talpha0(talpha0SEXP);
    Rcpp::traits::input_parameter< const arrayxd_external >::type tmu0(tmu0SEXP);
    Rcpp::traits::input_parameter< const arrayxd_external >::type tSiRiSr0(tSiRiSr0SEXP);
    Rcpp::traits::input_parameter< double >::type tolerance(toleranceSEXP);
    Rcpp::traits::input_parameter< int >::type itermax(itermaxSEXP);
    Rcpp::traits::input_parameter< Rcpp::LogicalVector >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< Rcpp::LogicalVector >::type lnz_tol(lnz_tolSEXP);
    rcpp_result_gen = Rcpp::wrap(rss_varbvsr_naive(SiRiS, sigma_beta, logodds, betahat, se, talpha0, tmu0, tSiRiSr0, tolerance, itermax, verbose, lnz_tol));
    return rcpp_result_gen;
END_RCPP
}
// rss_varbvsr_naive_sp
Rcpp::List rss_varbvsr_naive_sp(const sparseMatrix_external SiRiS, const double sigma_beta, const double logodds, const arrayxd_external betahat, const arrayxd_external se, const arrayxd_external talpha0, const arrayxd_external tmu0, const arrayxd_external tSiRiSr0, double tolerance, int itermax, Rcpp::LogicalVector verbose, Rcpp::LogicalVector lnz_tol);
RcppExport SEXP rssr_rss_varbvsr_naive_sp(SEXP SiRiSSEXP, SEXP sigma_betaSEXP, SEXP logoddsSEXP, SEXP betahatSEXP, SEXP seSEXP, SEXP talpha0SEXP, SEXP tmu0SEXP, SEXP tSiRiSr0SEXP, SEXP toleranceSEXP, SEXP itermaxSEXP, SEXP verboseSEXP, SEXP lnz_tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const sparseMatrix_external >::type SiRiS(SiRiSSEXP);
    Rcpp::traits::input_parameter< const double >::type sigma_beta(sigma_betaSEXP);
    Rcpp::traits::input_parameter< const double >::type logodds(logoddsSEXP);
    Rcpp::traits::input_parameter< const arrayxd_external >::type betahat(betahatSEXP);
    Rcpp::traits::input_parameter< const arrayxd_external >::type se(seSEXP);
    Rcpp::traits::input_parameter< const arrayxd_external >::type talpha0(talpha0SEXP);
    Rcpp::traits::input_parameter< const arrayxd_external >::type tmu0(tmu0SEXP);
    Rcpp::traits::input_parameter< const arrayxd_external >::type tSiRiSr0(tSiRiSr0SEXP);
    Rcpp::traits::input_parameter< double >::type tolerance(toleranceSEXP);
    Rcpp::traits::input_parameter< int >::type itermax(itermaxSEXP);
    Rcpp::traits::input_parameter< Rcpp::LogicalVector >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< Rcpp::LogicalVector >::type lnz_tol(lnz_tolSEXP);
    rcpp_result_gen = Rcpp::wrap(rss_varbvsr_naive_sp(SiRiS, sigma_beta, logodds, betahat, se, talpha0, tmu0, tSiRiSr0, tolerance, itermax, verbose, lnz_tol));
    return rcpp_result_gen;
END_RCPP
}
// rss_varbvsr_squarem
Rcpp::List rss_varbvsr_squarem(const Matrix_external SiRiS, const double sigma_beta, const double logodds, const arrayxd_external betahat, const arrayxd_external se, const arrayxd_external talpha0, const arrayxd_external tmu0, const arrayxd_external tSiRiSr0, double tolerance, int itermax, Rcpp::LogicalVector verbose, Rcpp::LogicalVector lnz_tol);
RcppExport SEXP rssr_rss_varbvsr_squarem(SEXP SiRiSSEXP, SEXP sigma_betaSEXP, SEXP logoddsSEXP, SEXP betahatSEXP, SEXP seSEXP, SEXP talpha0SEXP, SEXP tmu0SEXP, SEXP tSiRiSr0SEXP, SEXP toleranceSEXP, SEXP itermaxSEXP, SEXP verboseSEXP, SEXP lnz_tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Matrix_external >::type SiRiS(SiRiSSEXP);
    Rcpp::traits::input_parameter< const double >::type sigma_beta(sigma_betaSEXP);
    Rcpp::traits::input_parameter< const double >::type logodds(logoddsSEXP);
    Rcpp::traits::input_parameter< const arrayxd_external >::type betahat(betahatSEXP);
    Rcpp::traits::input_parameter< const arrayxd_external >::type se(seSEXP);
    Rcpp::traits::input_parameter< const arrayxd_external >::type talpha0(talpha0SEXP);
    Rcpp::traits::input_parameter< const arrayxd_external >::type tmu0(tmu0SEXP);
    Rcpp::traits::input_parameter< const arrayxd_external >::type tSiRiSr0(tSiRiSr0SEXP);
    Rcpp::traits::input_parameter< double >::type tolerance(toleranceSEXP);
    Rcpp::traits::input_parameter< int >::type itermax(itermaxSEXP);
    Rcpp::traits::input_parameter< Rcpp::LogicalVector >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< Rcpp::LogicalVector >::type lnz_tol(lnz_tolSEXP);
    rcpp_result_gen = Rcpp::wrap(rss_varbvsr_squarem(SiRiS, sigma_beta, logodds, betahat, se, talpha0, tmu0, tSiRiSr0, tolerance, itermax, verbose, lnz_tol));
    return rcpp_result_gen;
END_RCPP
}
// rss_varbvsr_squarem_sp
Rcpp::List rss_varbvsr_squarem_sp(const sparseMatrix_external SiRiS, const double sigma_beta, const double logodds, const arrayxd_external betahat, const arrayxd_external se, const arrayxd_external talpha0, const arrayxd_external tmu0, const arrayxd_external tSiRiSr0, double tolerance, int itermax, Rcpp::LogicalVector verbose, Rcpp::LogicalVector lnz_tol);
RcppExport SEXP rssr_rss_varbvsr_squarem_sp(SEXP SiRiSSEXP, SEXP sigma_betaSEXP, SEXP logoddsSEXP, SEXP betahatSEXP, SEXP seSEXP, SEXP talpha0SEXP, SEXP tmu0SEXP, SEXP tSiRiSr0SEXP, SEXP toleranceSEXP, SEXP itermaxSEXP, SEXP verboseSEXP, SEXP lnz_tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const sparseMatrix_external >::type SiRiS(SiRiSSEXP);
    Rcpp::traits::input_parameter< const double >::type sigma_beta(sigma_betaSEXP);
    Rcpp::traits::input_parameter< const double >::type logodds(logoddsSEXP);
    Rcpp::traits::input_parameter< const arrayxd_external >::type betahat(betahatSEXP);
    Rcpp::traits::input_parameter< const arrayxd_external >::type se(seSEXP);
    Rcpp::traits::input_parameter< const arrayxd_external >::type talpha0(talpha0SEXP);
    Rcpp::traits::input_parameter< const arrayxd_external >::type tmu0(tmu0SEXP);
    Rcpp::traits::input_parameter< const arrayxd_external >::type tSiRiSr0(tSiRiSr0SEXP);
    Rcpp::traits::input_parameter< double >::type tolerance(toleranceSEXP);
    Rcpp::traits::input_parameter< int >::type itermax(itermaxSEXP);
    Rcpp::traits::input_parameter< Rcpp::LogicalVector >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< Rcpp::LogicalVector >::type lnz_tol(lnz_tolSEXP);
    rcpp_result_gen = Rcpp::wrap(rss_varbvsr_squarem_sp(SiRiS, sigma_beta, logodds, betahat, se, talpha0, tmu0, tSiRiSr0, tolerance, itermax, verbose, lnz_tol));
    return rcpp_result_gen;
END_RCPP
}
// grid_search_rss_varbvsr_sp
Rcpp::DataFrame grid_search_rss_varbvsr_sp(const sparseMatrix_external SiRiS, const arrayxd_external sigma_beta, const arrayxd_external logodds, const arrayxd_external betahat, const arrayxd_external se, const arrayxd_external talpha0, const arrayxd_external tmu0, const arrayxd_external tSiRiSr0, double tolerance, int itermax, Rcpp::LogicalVector verbose, Rcpp::LogicalVector lnz_tol);
RcppExport SEXP rssr_grid_search_rss_varbvsr_sp(SEXP SiRiSSEXP, SEXP sigma_betaSEXP, SEXP logoddsSEXP, SEXP betahatSEXP, SEXP seSEXP, SEXP talpha0SEXP, SEXP tmu0SEXP, SEXP tSiRiSr0SEXP, SEXP toleranceSEXP, SEXP itermaxSEXP, SEXP verboseSEXP, SEXP lnz_tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const sparseMatrix_external >::type SiRiS(SiRiSSEXP);
    Rcpp::traits::input_parameter< const arrayxd_external >::type sigma_beta(sigma_betaSEXP);
    Rcpp::traits::input_parameter< const arrayxd_external >::type logodds(logoddsSEXP);
    Rcpp::traits::input_parameter< const arrayxd_external >::type betahat(betahatSEXP);
    Rcpp::traits::input_parameter< const arrayxd_external >::type se(seSEXP);
    Rcpp::traits::input_parameter< const arrayxd_external >::type talpha0(talpha0SEXP);
    Rcpp::traits::input_parameter< const arrayxd_external >::type tmu0(tmu0SEXP);
    Rcpp::traits::input_parameter< const arrayxd_external >::type tSiRiSr0(tSiRiSr0SEXP);
    Rcpp::traits::input_parameter< double >::type tolerance(toleranceSEXP);
    Rcpp::traits::input_parameter< int >::type itermax(itermaxSEXP);
    Rcpp::traits::input_parameter< Rcpp::LogicalVector >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< Rcpp::LogicalVector >::type lnz_tol(lnz_tolSEXP);
    rcpp_result_gen = Rcpp::wrap(grid_search_rss_varbvsr_sp(SiRiS, sigma_beta, logodds, betahat, se, talpha0, tmu0, tSiRiSr0, tolerance, itermax, verbose, lnz_tol));
    return rcpp_result_gen;
END_RCPP
}
// grid_search_rss_varbvsr
Rcpp::DataFrame grid_search_rss_varbvsr(const Matrix_external SiRiS, const arrayxd_external sigma_beta, const arrayxd_external logodds, const arrayxd_external betahat, const arrayxd_external se, const arrayxd_external talpha0, const arrayxd_external tmu0, const arrayxd_external tSiRiSr0, double tolerance, int itermax, Rcpp::LogicalVector verbose, Rcpp::LogicalVector lnz_tol);
RcppExport SEXP rssr_grid_search_rss_varbvsr(SEXP SiRiSSEXP, SEXP sigma_betaSEXP, SEXP logoddsSEXP, SEXP betahatSEXP, SEXP seSEXP, SEXP talpha0SEXP, SEXP tmu0SEXP, SEXP tSiRiSr0SEXP, SEXP toleranceSEXP, SEXP itermaxSEXP, SEXP verboseSEXP, SEXP lnz_tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Matrix_external >::type SiRiS(SiRiSSEXP);
    Rcpp::traits::input_parameter< const arrayxd_external >::type sigma_beta(sigma_betaSEXP);
    Rcpp::traits::input_parameter< const arrayxd_external >::type logodds(logoddsSEXP);
    Rcpp::traits::input_parameter< const arrayxd_external >::type betahat(betahatSEXP);
    Rcpp::traits::input_parameter< const arrayxd_external >::type se(seSEXP);
    Rcpp::traits::input_parameter< const arrayxd_external >::type talpha0(talpha0SEXP);
    Rcpp::traits::input_parameter< const arrayxd_external >::type tmu0(tmu0SEXP);
    Rcpp::traits::input_parameter< const arrayxd_external >::type tSiRiSr0(tSiRiSr0SEXP);
    Rcpp::traits::input_parameter< double >::type tolerance(toleranceSEXP);
    Rcpp::traits::input_parameter< int >::type itermax(itermaxSEXP);
    Rcpp::traits::input_parameter< Rcpp::LogicalVector >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< Rcpp::LogicalVector >::type lnz_tol(lnz_tolSEXP);
    rcpp_result_gen = Rcpp::wrap(grid_search_rss_varbvsr(SiRiS, sigma_beta, logodds, betahat, se, talpha0, tmu0, tSiRiSr0, tolerance, itermax, verbose, lnz_tol));
    return rcpp_result_gen;
END_RCPP
}
// wrap_rss_varbvsr_iter
Rcpp::List wrap_rss_varbvsr_iter(const Matrix_external SiRiS, const double sigma_beta, const double logodds, const arrayxd_external betahat, const arrayxd_external se, const arrayxd_external alpha, const arrayxd_external mu, const arrayxd_external SiRiSr, bool reverse);
RcppExport SEXP rssr_wrap_rss_varbvsr_iter(SEXP SiRiSSEXP, SEXP sigma_betaSEXP, SEXP logoddsSEXP, SEXP betahatSEXP, SEXP seSEXP, SEXP alphaSEXP, SEXP muSEXP, SEXP SiRiSrSEXP, SEXP reverseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Matrix_external >::type SiRiS(SiRiSSEXP);
    Rcpp::traits::input_parameter< const double >::type sigma_beta(sigma_betaSEXP);
    Rcpp::traits::input_parameter< const double >::type logodds(logoddsSEXP);
    Rcpp::traits::input_parameter< const arrayxd_external >::type betahat(betahatSEXP);
    Rcpp::traits::input_parameter< const arrayxd_external >::type se(seSEXP);
    Rcpp::traits::input_parameter< const arrayxd_external >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const arrayxd_external >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const arrayxd_external >::type SiRiSr(SiRiSrSEXP);
    Rcpp::traits::input_parameter< bool >::type reverse(reverseSEXP);
    rcpp_result_gen = Rcpp::wrap(wrap_rss_varbvsr_iter(SiRiS, sigma_beta, logodds, betahat, se, alpha, mu, SiRiSr, reverse));
    return rcpp_result_gen;
END_RCPP
}
// wrap_rss_varbvsr_iter_sp
Rcpp::List wrap_rss_varbvsr_iter_sp(const sparseMatrix_external SiRiS, const double sigma_beta, const double logodds, const arrayxd_external betahat, const arrayxd_external se, const arrayxd_external alpha, const arrayxd_external mu, const arrayxd_external SiRiSr, bool reverse);
RcppExport SEXP rssr_wrap_rss_varbvsr_iter_sp(SEXP SiRiSSEXP, SEXP sigma_betaSEXP, SEXP logoddsSEXP, SEXP betahatSEXP, SEXP seSEXP, SEXP alphaSEXP, SEXP muSEXP, SEXP SiRiSrSEXP, SEXP reverseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const sparseMatrix_external >::type SiRiS(SiRiSSEXP);
    Rcpp::traits::input_parameter< const double >::type sigma_beta(sigma_betaSEXP);
    Rcpp::traits::input_parameter< const double >::type logodds(logoddsSEXP);
    Rcpp::traits::input_parameter< const arrayxd_external >::type betahat(betahatSEXP);
    Rcpp::traits::input_parameter< const arrayxd_external >::type se(seSEXP);
    Rcpp::traits::input_parameter< const arrayxd_external >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const arrayxd_external >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const arrayxd_external >::type SiRiSr(SiRiSrSEXP);
    Rcpp::traits::input_parameter< bool >::type reverse(reverseSEXP);
    rcpp_result_gen = Rcpp::wrap(wrap_rss_varbvsr_iter_sp(SiRiS, sigma_beta, logodds, betahat, se, alpha, mu, SiRiSr, reverse));
    return rcpp_result_gen;
END_RCPP
}
// SiRSi
Eigen::SparseMatrix<double> SiRSi(const Eigen::MappedSparseMatrix<double>& R, const Eigen::VectorXd Si);
RcppExport SEXP rssr_SiRSi(SEXP RSEXP, SEXP SiSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MappedSparseMatrix<double>& >::type R(RSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type Si(SiSEXP);
    rcpp_result_gen = Rcpp::wrap(SiRSi(R, Si));
    return rcpp_result_gen;
END_RCPP
}
// SiRSi_d
Eigen::MatrixXd SiRSi_d(const Matrix_external R, const Eigen::VectorXd Si);
RcppExport SEXP rssr_SiRSi_d(SEXP RSEXP, SEXP SiSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Matrix_external >::type R(RSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type Si(SiSEXP);
    rcpp_result_gen = Rcpp::wrap(SiRSi_d(R, Si));
    return rcpp_result_gen;
END_RCPP
}
// genSymm
Eigen::SparseMatrix<double> genSymm(const Eigen::SparseMatrix<double>& R);
RcppExport SEXP rssr_genSymm(SEXP RSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::SparseMatrix<double>& >::type R(RSEXP);
    rcpp_result_gen = Rcpp::wrap(genSymm(R));
    return rcpp_result_gen;
END_RCPP
}
