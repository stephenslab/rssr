#ifndef INCLUDE_BVSR
#define INCLUDE_BVSR
#include <RcppEigen.h>





//typedef MSpMat::InnerIterator InIterMat;
//typedef Eigen::SparseVector<double> SpVec;
//typedef SpVec::InnerIterator InIterVec;


void rss_varbvsr_update (const double betahat, const double se, const double sigma_beta, const Eigen::ArrayXd &SiRiS_snp, 
                         Eigen::ArrayXd &SiRiSr, const double SiRiSr_snp, const double logodds, double &alpha, double &mu);

#endif
