
context("Alternate implementation")
library(Matrix)
# 
data(betahat)
# 
data(se)
# 
data(alpha_test)
data(mu_test)
data(R_shrink)
# 
betahat <- c(betahat)
se <- c(se)
mu_test <- c(mu_test)
alpha_test <- c(alpha_test)
t_SiRiS <- SiRSi(R_shrink,1/se)
SiRiS_f <- as.matrix(SiRSi(R_shrink,1/se))
p <- length(betahat)
SiRiSr=c(SiRiS_f%*%(alpha_test*mu_test))
SiRiS  <-as(SiRiS_f,"dgCMatrix")
# 
# 
sigb <- 1
logodds<--4

test_that("one iteration of alternate is equal to old implementation",{
  tres <- wrap_rss_varbvsr_iter_alt(t_SiRiS,sigb,logodds,betahat,se,alpha_test,mu_test,SiRiSr,F)
  mres <- wrap_rss_varbvsr_iter(t_SiRiS,sigb,logodds,betahat,se,alpha_test,mu_test,SiRiSr,F)
  
  expect_equal(tres,mres)
})

test_that("naive  alternate is equal to old implementation",{
  tres <- rss_varbvsr_naive(t_SiRiS,sigb,logodds,betahat,se,alpha_test,mu_test,SiRiSr,1e-4,100,T,F)
  mres <- rss_varbvsr_naive_alt(t_SiRiS,sigb,logodds,betahat,se,alpha_test,mu_test,SiRiSr,1e-4,100,T,F)
  
  expect_equal(tres,mres)
})
# mb <- microbenchmark(new=rss_varbvsr_naive_alt(t_SiRiS,sigb,logodds,betahat,se,alpha_test,mu_test,SiRiSr,1e-4,100,T,F),
#                      old=rss_varbvsr_naive(t_SiRiS,sigb,logodds,betahat,se,alpha_test,mu_test,SiRiSr,1e-4,100,T,F))
# 
# 
# test_that("mu and alpha are computed correctly",{
#   sigma_square <- ((se*se * (sigb*sigb)) / ((se*se) + (sigb*sigb)))[1]
#   mu_t <- wrap_compute_mu(betahat = betahat[1],se_square = (se*se)[1],sigma_square = sigma_square,alpha = alpha_test[1],mu = mu_test[1],SiRiSr_snp = SiRiSr[1])
#   expect_equal(mu_t,mres$mu1[1])
#   alpha_t <- wrap_compute_alpha(sigma_square =sigma_square,
#                                 sigma_beta=as.numeric(sigb),logodds = logodds,mu = mu_t,alpha = alpha_test[1] )
#   expect_equal(alpha_t,mres$alpha1[1])
# })


