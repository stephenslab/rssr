context("dense vs sparse (grid)")
library(Matrix)
library(purrr)
#Load the simulated data
data("betahat")
betahat <- c(betahat)
data("se")
se <- c(se)
data("alpha_test")
data("mu_test")
data("R_shrink")
#Convert the rowvector to a column vector
mu_test <- t(t(mu_test))
alpha_test <- t(t(alpha_test))
# R_panel <-as.matrix(R_panel)

#Generate SiRiS as both dense and sparse matrices 
SiRiS <- SiRSi(R_shrink,1/se)
SiRiS_f <- SiRSi_d(as.matrix(R_shrink),1/se)
#SiRiS_f <- as.matrix(SiRSi(R_shrink,1/se))
#SiRiS <-as(SiRiS_f,"dgCMatrix")
p <- length(betahat)
SiRiSr=c(SiRiS_f%*%(alpha_test*mu_test))


test_that("grid implementations are identical",{
  sigb <- c(1,1.1,1.4)
  logodds <- c(-4.3,-4.4,-4.5)
  my_grid_sp <- grid_search_rss_varbvsr_sparse(SiRiS = SiRiS,sigma_beta = sigb,logodds = logodds,betahat = betahat,se = se,alpha0 =  alpha_test,mu0 =  mu_test,SiRiSr0 = SiRiSr,tolerance = 1e-4,itermax = 100,lnz_tol = F)
  #my_grid_d <- grid_search_rss_varbvsr_dense(SiRiS = SiRiS_f,sigma_beta = sigb,logodds = -4.6,betahat = betahat,se = se,talpha0 =  alpha_test,tmu0 =   mu_test,tSiRiSr0 = SiRiSr,tolerance = 1e-4,itermax = 100,lnz_tol = F,verbose=T)
  my_grid_dd <- grid_search_rss_varbvsr_dense(SiRiS = SiRiS_f,sigma_beta = sigb,logodds = logodds,betahat = betahat,se = se,alpha0 =  alpha_test,mu0 =  mu_test,SiRiSr0 = SiRiSr,tolerance = 1e-4,itermax = 100,lnz_tol = F)

  # mb <- microbenchmark::microbenchmark(sparse=grid_search_rss_varbvsr_sparse(SiRiS = SiRiS,sigma_beta = sigb,logodds = logodds,betahat = betahat,se = se,alpha0 =  alpha_test,mu0 =  mu_test,SiRiSr0 = SiRiSr,tolerance = 1e-4,itermax = 100,lnz_tol = F),
                 # dense= grid_search_rss_varbvsr_dense(SiRiS = SiRiS_f,sigma_beta = sigb,logodds = logodds,betahat = betahat,se = se,alpha0 =  alpha_test,mu0 =  mu_test,SiRiSr0 = SiRiSr,tolerance = 1e-4,itermax = 100,lnz_tol = F))
  expect_equivalent(my_grid_dd,my_grid_sp)
})



