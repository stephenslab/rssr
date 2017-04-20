context("array vs not-array")

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
SiRiS_f <- SiRSi_d(as.matrix(R_shrink),1/se)
#SiRiS_f <- as.matrix(SiRSi(R_shrink,1/se))
SiRiS <-as(SiRiS_f,"dgCMatrix")
p <- length(betahat)
SiRiSr=c(SiRiS_f%*%(alpha_test*mu_test))


test_that("Naive implementations are identical",{
  sigb <- 1
  
  my_naive_a <- rss_varbvsr_naive_array(SiRiS = SiRiS_f,sigma_beta = sigb,logodds = -4.6,betahat = betahat,se = se,talpha0 =  alpha_test,tmu0 =  mu_test,tSiRiSr0 = SiRiSr,tolerance = 1e-4,itermax = 100,verbose = T,lnz_tol = F)
  my_naive_s <- rss_varbvsr_naive(SiRiS = SiRiS_f,sigma_beta = sigb,logodds = -4.6,betahat = betahat,se = se,talpha0 =  alpha_test,tmu0 =  mu_test,tSiRiSr0 = SiRiSr,tolerance = 1e-4,itermax = 100,verbose = T,lnz_tol = F)
  expect_equivalent(c(my_naive_a$lnZ),c(my_naive_s$lnZ))
  expect_equivalent(c(my_naive_a$mu),c(my_naive_s$mu))
  expect_equivalent(c(my_naive_a$alpha),c(my_naive_s$alpha))
  expect_equivalent(my_naive_a$iter,my_naive_s$iter)
  expect_equivalent(my_naive_a$maxerr,my_naive_s$maxerr)
})


test_that("SQUAREM updates are identical",{
  sigb <- 1  
  logodds <- -3
  
  my_results_s <- wrap_rss_varbvsr_iter(SiRiS = SiRiS_f,
                                      sigma_beta=sigb,
                                      logodds=logodds,
                                      betahat = betahat,
                                      se = se,alpha = alpha_test,mu = mu_test,SiRiSr = SiRiSr,reverse = F)
  my_results_a <- wrap_rss_varbvsr_iter_array(SiRiS = SiRiS_f,
                            sigma_beta=sigb,
                            logodds=logodds,
                            betahat = betahat,
                            se = se,alpha = alpha_test,mu = mu_test,SiRiSr = SiRiSr,reverse = F)
  
  
  expect_equivalent(my_results_s$alpha1,c(my_results_a$alpha1))
  expect_equivalent(my_results_s$mu1,c(my_results_a$mu1))
  expect_equivalent(my_results_s$SiRiSr,c(my_results_a$SiRiSr))
})


test_that("grid optimization over logodds works as expected",{
  sigb <- 1
  log10oddsvec <- seq(-6,-1,0.5)
  logoddsvec <- log10oddsvec*log(10)

  
  mr_a <- grid_search_rss_varbvsr_array(SiRiS=SiRiS_f,sigma_beta =1,logodds=logoddsvec,betahat=betahat,se=se,talpha0=alpha_test,tmu0=mu_test,tSiRiSr0=SiRiSr,1e-4,100,F,F)  
  mr_s <- grid_search_rss_varbvsr(      SiRiS=SiRiS_f,sigma_beta =1,logodds=logoddsvec,betahat=betahat,se=se,talpha0=alpha_test,tmu0=mu_test,tSiRiSr0=SiRiSr,1e-4,100,F,F)  
  expect_equal(mr_a,mr_s)
})


test_that("2d grid optimization over sigb and logodds works the same",{
  log10oddsvec <- seq(-3.1,-2.1,length.out = 5)
  logoddsvec <- log10oddsvec*log(10)
  sigb <- seq(0.8,1.2,length.out = 4)
  mr_grid_s <- grid_search_rss_varbvsr(talpha0=alpha_test,
                                           tmu0=mu_test,betahat=betahat,
                                           se=se,
                                           SiRiS=SiRiS_f,
                                           sigma_beta =sigb,
                                           logodds=logoddsvec,
                                           verbose=F,
                                           tSiRiSr0=SiRiSr,itermax=100,tolerance=1e-4,lnz_tol=F)
  mr_grid_a <- grid_search_rss_varbvsr_array(talpha0=alpha_test,
                                           tmu0=mu_test,betahat=betahat,
                                           se=se,
                                           SiRiS=SiRiS_f,
                                           sigma_beta =sigb,
                                           logodds=logoddsvec,
                                           verbose=F,
                                           tSiRiSr0=SiRiSr,itermax=100,tolerance=1e-4,lnz_tol=F)
  
  
  expect_equal(mr_grid_a,mr_grid_s)
  }
)



