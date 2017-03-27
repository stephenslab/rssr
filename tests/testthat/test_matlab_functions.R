context("MATLAB_scripts")
library(RcppOctave)
#library(rssr)
library(Matrix)
library(testthat)
#Find the location of the .m files 
mfile <- system.file("m_files/run_install.m",package="rssr")
mdir <- system.file("m_files",package="rssr")

#change to the directory with the .m files in Octave
.CallOctave('cd',mdir)
o_source("run_install.m")

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
t_SiRiS <- SiRSi(R_shrink,1/se)
SiRiS_f <- as.matrix(SiRSi(R_shrink,1/se))
SiRiS <-as(SiRiS_f,"dgCMatrix")
p <- length(betahat)
SiRiSr=c(SiRiS_f%*%(alpha_test*mu_test))





#Call with Octave and RSSR



test_that("Single RSS update of alpha,mu and SiRiSr are approximately equal",{
  I <- 1:p
  rI <- p:1
  sigb <- 1
  logodds <- -3
  res <- .CallOctave('wrap_rss_varbvsr_update',SiRiS_f,sigb,logodds,betahat,se,alpha_test,mu_test,SiRiSr,I)
  mres <- wrap_rss_varbvsr_iter(t_SiRiS,sigb,logodds,betahat,se,alpha_test,mu_test,SiRiSr,F)
  expect_equal(c(res$alpha1),c(mres$alpha1),tolerance=1e-8)
  expect_equal(c(res$mu1),c(mres$mu1),tolerance=1e-8)
  expect_equal(c(res$SiRiSr),c(mres$SiRiSr),tolerance=1e-8)})


test_that("Single RSS update is the same when computed backwards",{
  I <- 1:p
  rI <- p:1
  sigb <- 1
  logodds <- -3
  rres <- .CallOctave('wrap_rss_varbvsr_update',SiRiS_f,sigb,logodds,betahat,se,alpha_test,mu_test,SiRiSr,rI)
  rmres <- wrap_rss_varbvsr_iter(t_SiRiS,sigb,logodds,betahat,se,alpha_test,mu_test,SiRiSr,T)
  expect_equal(c(rres$alpha1),c(rmres$alpha1),tolerance=1e-8)
  expect_equal(c(rres$mu1),c(rmres$mu1),tolerance=1e-8)
  expect_equal(c(rres$SiRiSr),c(rmres$SiRiSr),tolerance=1e-8)})



test_that("gamma integral is calcualted correctly",{
          logodds <- -3
          expect_equal(intgamma(logodds,alpha_test),
                       .CallOctave('intgamma',logodds,alpha_test))})

test_that("integral of variational lower bound is computed correctly",{
  sigb <- 1
  s_test=t(t(se*se*(sigb*sigb)/(se*se+sigb)))
  expect_equal(intklbeta_rssbvsr(alpha_test,mu_test,s_test,sigb),
               .CallOctave('intklbeta_rssbvsr',t(t(alpha_test)),t(t(mu_test)),s_test,sigb))
          })

test_that("betavar works the same",{
  sigb <- 1
  s_test=t(t(se*se*(sigb*sigb)/(se*se+sigb)))
  expect_equal(betavar(alpha_test,mu_test,s_test),
               c(.CallOctave('betavar',alpha_test,mu_test,s_test)))
  })


test_that("SiRiS is generated equivalently",{
  t_SiRiS = .CallOctave('gen_SiRiS',as.matrix(R_shrink),se)
  m_SiRiS <- as.matrix(SiRSi(R_shrink,1/se))
  attr(m_SiRiS,"dimnames") <- NULL
  expect_equivalent(t_SiRiS,m_SiRiS)
})



test_that("Naive implementations are identical",{
  sigb <- 1

  naive_results <- .CallOctave('wrap_rss_varbvsr_naive',t(t(betahat)),t(t(se)),SiRiS_f,sigb,-4.6,t(alpha_test),t(mu_test))
  my_naive <- rss_varbvsr_naive(SiRiS = SiRiS,sigma_beta = sigb,logodds = -4.6,betahat = betahat,se = se,talpha0 =  alpha_test,tmu0 =  mu_test,tSiRiSr0 = SiRiSr,tolerance = 1e-4,itermax = 100,verbose = T,lnz_tol = F)
  expect_equivalent(naive_results$lnZ,my_naive$lnZ)
  expect_equivalent(c(naive_results$mu),c(my_naive$mu))
  expect_equivalent(c(naive_results$alpha),c(my_naive$alpha))
  expect_equivalent(naive_results$info$iter,my_naive$iter)
  expect_equivalent(naive_results$info$max_err,my_naive$maxerr)
})



test_that("SQUAREM updates are identical",{
  sigb <- 1  
  logodds <- -3
  mat_results <- .CallOctave('wrap_rss_varbvsr_squarem',t(t(betahat)),t(t(se)),SiRiS_f,sigb,logodds,t(alpha_test),t(mu_test),1e-4)
  my_results <- rss_varbvsr_squarem(SiRiS = SiRiS,
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
  
  expect_equivalent(my_results$lnZ,mat_results$lnZ)
  expect_equivalent(my_results$mu,c(mat_results$mu))
  expect_equivalent(my_results$alpha,c(mat_results$alpha))
  expect_equivalent(my_results$iter,mat_results$info$iter)
  expect_equivalent(my_results$max_err,mat_results$info$maxerr)
})



test_that("grid optimization over logodds works as expected",{
  sigb <- 1
  log10oddsvec <- seq(-6,-1,0.5)
  logoddsvec <- log10oddsvec*log(10)
  pm <- .CallOctave('grid_rss_varbvsr_logodds',t(t(betahat)),t(t(se)),SiRiS_f,sigb,log10oddsvec,t(alpha_test),t(mu_test))
  mr <- grid_search_rss_varbvsr(SiRiS=SiRiS,sigma_beta =1,logodds=logoddsvec,betahat=betahat,se=se,talpha0=alpha_test,tmu0=mu_test,tSiRiSr0=SiRiSr,1e-4,100,F,F)  
  pi_mean <- marg_pi(log10odds = log10oddsvec,c(mr$lnZ))  
  expect_equal(c(pi_mean),pm$pi_mean)
  })

test_that("2d grid optimization over sigb and logodds works as in MATLAB",{
  log10oddsvec <- seq(-3.1,-2.1,length.out = 5)
  logoddsvec <- log10oddsvec*log(10)
  sigb <- seq(0.8,1.2,length.out = 5)
  pm <- .CallOctave('grid_rssr_varbvsr',t(t(betahat)),t(t(se)),SiRiS_f,sigb,log10oddsvec,t(alpha_test),t(mu_test))
  mr_grid <- grid_search_rss_varbvsr(talpha0=alpha_test,
                                     tmu0=mu_test,betahat=betahat,
                                     se=se,
                                     SiRiS=SiRiS,
                                     sigma_beta =sigb,
                                     logodds=logoddsvec,
                                     verbose=F,
                                     tSiRiSr0=SiRiSr,itermax=100,tolerance=1e-4,lnz_tol=F)
  expect_equal(c(mr_grid$lnZ),c(pm))})






# hfile <- "/media/nwknoblauch/Data/GTEx/1kg_LD/EUR.chr1_1_1kg.h5"
# SiRiS_f <- as.matrix(gen_SiRSi(hfile))
# #R <- as.matrix(read_ccs_h5("/media/nwknoblauch/Data/GTEx/1kg_LD/EUR.chr1_1_1kg.h5",groupname = "R"))
# betahat <- read_vec(hfile,"betahat")
# se <- read_vec(hfile,"se")
# p <- length(se)
# alpha <- ralpha(p)
# mu <- rmu(p)
# sigb <- 0.058
# logodds <- -2.9/log(10)
# mat_results <- .CallOctave('wrap_rss_varbvsr_squarem',t(t(betahat)),t(t(se)),SiRiS_f,sigb,logodds,t(alpha),t(mu))



