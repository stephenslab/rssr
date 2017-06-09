context("MATLAB_scripts")
# library(rssr)
library(Matrix)
library(rssr)
#Load the simulated data
data("betahat")
betahat <- c(betahat)
data("se")
se <- c(se)
data("alpha_test")
data("mu_test")
data("R_shrink")
# evd <- eigen(R_shrink)
#Convert the rowvector to a column vector
mu_test <- t(t(mu_test))
alpha_test <- t(t(alpha_test))
Sample_Size <- 1458
# R_panel <-as.matrix(R_panel)

#Generate SiRiS as both dense and sparse matrices 
SiRiS <- SiRSi(R_shrink,1/se)
SiRiS_f <- SiRSi_d(as.matrix(R_shrink),1/se)
#SiRiS_f <- as.matrix(SiRSi(R_shrink,1/se))
#SiRiS <-as(SiRiS_f,"dgCMatrix")
p <- length(betahat)
SiRiSr=c(SiRiS_f%*%(alpha_test*mu_test))







test_that("SQUAREM step size adjustment works",{
  
  sigb <- 1
  logodds <- -3
  
  data("matlab_update")
  mat_results <- matlab_update
  t_adjust <-  wrap_squarem_adjust_prep(SiRiS = SiRiS_f,
                                        sigma_beta=sigb,
                                        logodds=logodds,
                                        betahat = betahat,
                                        se = se,
                                        talpha = alpha_test,
                                        tmu = mu_test,
                                        tSiRiSr = SiRiSr,
                                        tolerance = 1e-4,
                                        itermax=0,
                                        lnz_tol=F)
  
  expect_equal(t_adjust$mtp,mat_results$mtp)
  expect_equal(t_adjust$alpha0,c(mat_results$alpha0))
  expect_equal(t_adjust$alpha1,c(mat_results$alpha1))
  expect_equal(t_adjust$alpha2,c(mat_results$alpha2))
  expect_equal(t_adjust$alpha,c(mat_results$alpha))
  expect_equal(t_adjust$mu,c(mat_results$mu))
  expect_equal(t_adjust$mu0,c(mat_results$mu0))
  expect_equal(t_adjust$mu,c(mat_results$mu))
})



test_that("Single RSS update of alpha,mu and SiRiSr are approximately equal",{


  I <- 1:p
  rI <- p:1
  sigb <- 1
  logodds <- -3
  data("matlab_varbvsr_update_1")
  res <- matlab_varbvsr_update_1
  mres <- wrap_rss_varbvsr_iter_sp(SiRiS,sigb,logodds,betahat,se,alpha_test,mu_test,SiRiSr,F)
  expect_equal(c(res$alpha1),c(mres$alpha1),tolerance=1e-8)
  expect_equal(c(res$mu1),c(mres$mu1),tolerance=1e-8)
  expect_equal(c(res$SiRiSr),c(mres$SiRiSr),tolerance=1e-8)})


test_that("Single RSS update is the same when computed backwards",{
  I <- 1:p
  rI <- p:1
  sigb <- 1
  logodds <- -3
  data("matlab_varbvsr_update_2")
  rres <- matlab_varbvsr_update_2
  rmres <- wrap_rss_varbvsr_iter_sp(SiRiS,sigb,logodds,betahat,se,alpha_test,mu_test,SiRiSr,T)
  expect_equal(c(rres$alpha1),c(rmres$alpha1),tolerance=1e-8)
  expect_equal(c(rres$mu1),c(rmres$mu1),tolerance=1e-8)
  expect_equal(c(rres$SiRiSr),c(rmres$SiRiSr),tolerance=1e-8)})


test_that("SiRiS is generated equivalently",{
  data("matlab_SiRiS")
t_SiRiS <- matlab_SiRiS
  m_SiRiS <- as.matrix(SiRSi(R_shrink,1/se))
  attr(m_SiRiS,"dimnames") <- NULL
  expect_equivalent(t_SiRiS,m_SiRiS)
})






test_that("Naive implementations are identical",{
  sigb <- 1
  data("matlab_varbvsr_naive_1")
  naive_results <- matlab_varbvsr_naive_1
  my_naive <- rss_varbvsr_naive_sp(SiRiS = SiRiS,sigma_beta = sigb,logodds = -4.6,betahat = betahat,se = se,talpha0 =  alpha_test,tmu0 =  mu_test,tSiRiSr0 = SiRiSr,tolerance = 1e-4,itermax = 100,verbose = T,lnz_tol = F)
  expect_equivalent(naive_results$lnZ,my_naive$lnZ)
  expect_equivalent(c(naive_results$mu),c(my_naive$mu))
  expect_equivalent(c(naive_results$alpha),c(my_naive$alpha))
  expect_equivalent(naive_results$info$iter,my_naive$iter)
  expect_equivalent(naive_results$info$max_err,my_naive$maxerr)
})





test_that("SQUAREM updates are identical",{
  sigb <- 1  
  logodds <- -3
  data("matlab_varbvsr_squarem_1")
  mat_results_2 <- matlab_varbvsr_squarem_1
  # mat_results_1up <- matlab_varbvsr_squarem_1_up
  my_results_2 <- rss_varbvsr_squarem_sp(SiRiS = SiRiS,
                                         sigma_beta=sigb,
                                         logodds=logodds,
                                         betahat = betahat,
                                         se = se,
                                         talpha0 = alpha_test,
                                         tmu0 = mu_test,
                                         tSiRiSr0 = SiRiSr,
                                         tolerance = 1e-4,
                                         itermax=100,
                                         verbose=T,
                                         lnz_tol=F)
  expect_equal(my_results_2$alpha,c(mat_results_2$alpha))
  expect_equal(my_results_2$mu,c(mat_results_2$mu))

  
  
  my_results <- rss_varbvsr_squarem_sp(SiRiS = SiRiS,
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

  
  expect_equivalent(my_results$lnZ,mat_results_2$lnZ)
  expect_equivalent(my_results$mu,c(mat_results_2$mu))
  expect_equivalent(my_results$alpha,c(mat_results_2$alpha))
  expect_equivalent(my_results$iter,mat_results_2$info$iter)
  expect_equivalent(my_results$max_err,mat_results_2$info$maxerr)
})



test_that("grid optimization over logodds works as expected",{
  sigb <- 1
  log10oddsvec <- seq(-6,-1,0.5)
  logoddsvec <- log10oddsvec*log(10)
  data("matlab_grid_logodds")
  pm <- matlab_grid_logodds
  paramdf <- list(sigb=sigb,logodds=logoddsvec) %>% purrr::cross_d()
  mr <- grid_search_rss_varbvsr_sp(SiRiS=SiRiS,sigma_beta =paramdf$sigb,logodds=paramdf$logodds,betahat=betahat,se=se,talpha0=alpha_test,tmu0=mu_test,tSiRiSr0=SiRiSr,1e-4,100,F,F)  

  mr_w <- normalizeLogWeights(mr$lnZ)
  pi <- exp(mr$logodds)/(1+exp(mr$logodds))
  R_pi_mean <- sum(pi*mr_w)
  expect_equal(R_pi_mean,pm$pi_mean)
})

test_that("2d grid optimization over sigb and logodds works as in MATLAB",{
  log10oddsvec <- seq(-3.1,-2.1,length.out = 5)
  logoddsvec <- log10oddsvec*log(10)
  sigb <- seq(0.8,1.2,length.out = 5)
  
  
  data("matlab_grid_logodds_sigb")
  pm <- matlab_grid_logodds_sigb
  
  paramdf <- list(sigb=sigb,logodds=logoddsvec) %>% purrr::cross_d()
  mr_grid <- grid_search_rss_varbvsr_sp(talpha0=alpha_test,
                                        tmu0=mu_test,betahat=betahat,
                                        se=se,
                                        SiRiS=SiRiS,
                                        sigma_beta =paramdf$sigb,
                                        logodds=paramdf$logodds,
                                        verbose=F,
                                        tSiRiSr0=SiRiSr,itermax=100,tolerance=1e-4,lnz_tol=F)
  mr_grid$pve <- mr_grid$pve/Sample_Size
  expect_equal(c(mr_grid$lnZ),c(t(pm)))}
)





test_that("2d parallel grid optimization over sigb and logodds works as in old version",{
  library(dplyr)
  log10oddsvec <- seq(-3.1,-2.1,length.out = 10)
  logoddsvec <- log10oddsvec*log(10)
  sigb <- seq(0.8,1.2,length.out = 10)
  
  
  # data("matlab_grid_logodds_sigb")
  # pm <- matlab_grid_logodds_sigb
  paramdf <- list(sigb=sigb,logodds=logoddsvec) %>% purrr::cross_d()
  SiRiSf <- as.matrix(SiRiS)
  SiRiS_f <- SiRSi_d(as.matrix(R_shrink),1/se)
  SiRiSr <- c(SiRiS_f%*%(alpha_test*mu_test))
  alpha_test <- c(alpha_test)
  mu_test <- c(mu_test)
  mr_grid_d <- grid_search_rss_varbvsr_alt(talpha0=alpha_test,
                                           tmu0=mu_test,
                                           betahat=betahat,
                                           se=se,
                                           SiRiS=SiRiS_f,
                                           sigma_beta =paramdf$sigb,
                                           logodds=paramdf$logodds,
                                           tSiRiSr0=SiRiSr,
                                           itermax=100L,
                                           tolerance=1e-4,
                                           lnz_tol=F,
                                           n=1L)
  
  mr_grid_o <- grid_search_rss_varbvsr(talpha0=alpha_test,
                                       tmu0=mu_test,
                                       betahat=betahat,
                                       se=se,
                                       SiRiS=SiRiS_f,
                                       sigma_beta =paramdf$sigb,
                                       logodds=paramdf$logodds,
                                       tSiRiSr0=SiRiSr,
                                       itermax=100L,
                                       tolerance=1e-4,
                                       lnz_tol=F,verbose = F)
  testthat::expect_equal(mr_grid_d,mr_grid_o)
  
  
  fmbl <- list()  
  pb <- progress::progress_bar$new(total=50)
  for(i in 1:50){
    
    fmbl[[i]] <- data.frame(microbenchmark::microbenchmark(new_tbb=grid_search_rss_varbvsr_alt(talpha0=alpha_test,
                                                                                               tmu0=mu_test,
                                                                                               betahat=betahat,
                                                                                               se=se,
                                                                                               SiRiS=SiRiS_f,
                                                                                               sigma_beta =paramdf$sigb,
                                                                                               logodds=paramdf$logodds,
                                                                                               tSiRiSr0=SiRiSr,
                                                                                               itermax=100L,
                                                                                               tolerance=1e-4,
                                                                                               lnz_tol=F,
                                                                                               n=1L,grainsize = i))) %>% mutate(grainsize=i)
    pb$tick()
  }
  ores <- microbenchmark(old_tbb=grid_search_rss_varbvsr(talpha0=alpha_test,
                                                                    tmu0=mu_test,
                                                                    betahat=betahat,
                                                                    se=se,
                                                                    SiRiS=SiRiS_f,
                                                                    sigma_beta =paramdf$sigb,
                                                                    logodds=paramdf$logodds,
                                                                    tSiRiSr0=SiRiSr,
                                                                    itermax=100L,
                                                                    tolerance=1e-4,
                                                                    lnz_tol=F,verbose = F),
                                    new_tbb=grid_search_rss_varbvsr_alt(talpha0=alpha_test,
                                                                        tmu0=mu_test,
                                                                        betahat=betahat,
                                                                        se=se,
                                                                        SiRiS=SiRiS_f,
                                                                        sigma_beta =paramdf$sigb,
                                                                        logodds=paramdf$logodds,
                                                                        tSiRiSr0=SiRiSr,
                                                                        itermax=100L,
                                                                        tolerance=1e-4,
                                                                        lnz_tol=F,
                                                                        n=1L,grainsize = 1))

  afmb <- bind_rows(fmbl)
  head()
  ggplot(afmb,aes(x=grainsize,y=time,group=grainsize))+geom_boxplot()
  
  expect_equal(mr_grid_d,mr_grid_o)
  
  
  
}
)









