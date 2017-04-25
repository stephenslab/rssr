context("dense vs sparse")

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


test_that("Naive implementations are identical",{
  sigb <- 1
  
  my_naive_sp <- rss_varbvsr_naive_sp(SiRiS = SiRiS,sigma_beta = sigb,logodds = -4.6,betahat = betahat,se = se,talpha0 =  alpha_test,tmu0 =  mu_test,tSiRiSr0 = SiRiSr,tolerance = 1e-4,itermax = 100,verbose = T,lnz_tol = F)
  my_naive_d <- rss_varbvsr_naive(SiRiS = SiRiS_f,sigma_beta = sigb,logodds = -4.6,betahat = betahat,se = se,talpha0 =  alpha_test,tmu0 =  mu_test,tSiRiSr0 = SiRiSr,tolerance = 1e-4,itermax = 100,verbose = T,lnz_tol = F)
  expect_equivalent(my_naive_d$lnZ,my_naive_sp$lnZ)
  expect_equivalent(c(my_naive_d$mu),c(my_naive_sp$mu))
  expect_equivalent(c(my_naive_d$alpha),c(my_naive_sp$alpha))
  expect_equivalent(my_naive_d$iter,my_naive_sp$iter)
  expect_equivalent(my_naive_d$maxerr,my_naive_sp$maxerr)
})


test_that("SQUAREM updates are identical",{
  sigb <- 1  
  logodds <- -3
  my_results_sp <- rss_varbvsr_squarem_sp(SiRiS = SiRiS,
                                          sigma_beta=sigb,
                                          logodds=logodds,
                                          betahat = betahat,
                                          se = se,
                                          talpha0 = alpha_test,
                                          tmu0 = mu_test,
                                          tSiRiSr0 = SiRiSr,
                                          tolerance = 1e-4,
                                          itermax=900,
                                          verbose=T,
                                          lnz_tol=F)
  my_results_d <- rss_varbvsr_squarem(SiRiS = SiRiS_f,
                                          sigma_beta=sigb,
                                          logodds=logodds,
                                          betahat = betahat,
                                          se = se,
                                          talpha0 = alpha_test,
                                          tmu0 = mu_test,
                                          tSiRiSr0 = SiRiSr,
                                          tolerance = 1e-4,
                                          itermax=900,
                                          verbose=T,
                                          lnz_tol=F)
  
  
  expect_equivalent(my_results_d,my_results_sp)
})


test_that("grid optimization over logodds works as expected",{
  sigb <- 1
  log10oddsvec <- seq(-6,-1,0.5)
  logoddsvec <- log10oddsvec*log(10)
  mr_s <- grid_search_rss_varbvsr_sp(SiRiS=SiRiS,sigma_beta =1,logodds=logoddsvec,betahat=betahat,se=se,talpha0=alpha_test,tmu0=mu_test,tSiRiSr0=SiRiSr,1e-4,100,F,F)  
  mr_d <- grid_search_rss_varbvsr(SiRiS=SiRiS_f,sigma_beta =1,logodds=logoddsvec,betahat=betahat,se=se,talpha0=alpha_test,tmu0=mu_test,tSiRiSr0=SiRiSr,1e-4,100,F,F)  
  expect_equal(mr_s$logodds,mr_d$logodds)
  expect_equal(mr_s$sigb,mr_d$sigb)
  expect_equal(mr_s,mr_d)
})



test_that("dense is faster than sparse",{
  sigb <- 1
  log10oddsvec <- seq(-6,-1,0.5)
  logoddsvec <- log10oddsvec*log(10)
  time_s <- system.time(mr_s <- grid_search_rss_varbvsr_sp(SiRiS=SiRiS,sigma_beta =1,logodds=logoddsvec,betahat=betahat,se=se,talpha0=alpha_test,tmu0=mu_test,tSiRiSr0=SiRiSr,1e-4,100,F,F)  )
  time_d <- system.time(mr_d <- grid_search_rss_varbvsr(SiRiS=SiRiS_f,sigma_beta =1,logodds=logoddsvec,betahat=betahat,se=se,talpha0=alpha_test,tmu0=mu_test,tSiRiSr0=SiRiSr,1e-4,100,F,F)  )
  expect_lt(time_d["elapsed"],time_s["elapsed"])
})

test_that("2d grid optimization over sigb and logodds works the same",{
  log10oddsvec <- seq(-3.1,-2.1,length.out = 5)
  logoddsvec <- log10oddsvec*log(10)
  sigb <- seq(0.8,1.2,length.out = 5)
  mr_grid_sp <- grid_search_rss_varbvsr_sp(talpha0=alpha_test,
                                           tmu0=mu_test,betahat=betahat,
                                           se=se,
                                           SiRiS=SiRiS,
                                           sigma_beta =sigb,
                                           logodds=logoddsvec,
                                           verbose=F,
                                           tSiRiSr0=SiRiSr,itermax=100,tolerance=1e-4,lnz_tol=F)
  mr_grid_d <- grid_search_rss_varbvsr(talpha0=alpha_test,
                                           tmu0=mu_test,betahat=betahat,
                                           se=se,
                                           SiRiS=SiRiS_f,
                                           sigma_beta =sigb,
                                           logodds=logoddsvec,
                                           verbose=F,
                                           tSiRiSr0=SiRiSr,itermax=100,tolerance=1e-4,lnz_tol=F)
  
  
  expect_equal(c(mr_grid_d),c(mr_grid_sp))}
)



