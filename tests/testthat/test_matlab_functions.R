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

sigb <- 1
logodds <- -3
s_test=t(t(se*se*(sigb*sigb)/(se*se+sigb)))
I <- 1:p

#Call with Octave and RSSR
res <- .CallOctave('wrap_rss_varbvsr_update',SiRiS_f,sigb,logodds,betahat,se,alpha_test,mu_test,SiRiSr,p)
mres <- wrap_rss_varbvsr_iter(t_SiRiS,sigb,logodds,betahat,se,alpha_test,mu_test,SiRiSr,F)


test_that("Single RSS update of alpha,mu and SiRiSr are approximately equal",{
  expect_equal(c(res$alpha1),c(mres$alpha1),tolerance=1e-8)
  expect_equal(c(res$mu1),c(mres$mu1),tolerance=1e-8)
  expect_equal(c(res$SiRiSr),c(mres$SiRiSr),tolerance=1e-8)})


test_that("gamma integral is calcualted correctly",
          expect_equal(intgamma(logodds,alpha_test),
                       .CallOctave('intgamma',logodds,alpha_test)))

test_that("integral of variational lower bound is computed correctly",
          expect_equal(intklbeta_rssbvsr(alpha_test,mu_test,s_test,sigb),
                       .CallOctave('intklbeta_rssbvsr',alpha_test,mu_test,s_test,sigb)))

test_that("betavar works the same",
          expect_equal(betavar(alpha_test,mu_test,s_test),
                       c(.CallOctave('betavar',alpha_test,mu_test,s_test))))

t_SiRiS = .CallOctave('gen_SiRiS',as.matrix(R_shrink),se)
m_SiRiS <- as.matrix(SiRSi(R_shrink,1/se))
attr(m_SiRiS,"dimnames") <- NULL
test_that("SiRiS is generated equivalently",expect_equivalent(t_SiRiS,m_SiRiS))
rm(m_SiRiS,t_SiRiS)

mat_results <- .CallOctave('wrap_rss_varbvsr_squarem',t(t(betahat)),t(t(se)),SiRiS_f,sigb,logodds,t(alpha_test),t(mu_test))
naive_results <- .CallOctave('wrap_rss_varbvsr_naive',t(t(betahat)),t(t(se)),SiRiS_f,sigb,-4.6,t(alpha_test),t(mu_test))
my_naive <- rss_varbvsr_naive(SiRiS = SiRiS,sigma_beta = sigb,logodds = -4.6,betahat = betahat,se = se,alpha0 = alpha_test,mu0 = mu_test,SiRiSr0 = SiRiSr,tolerance = 1e-4)


test_that("Naive implementations are identical",{
  expect_equivalent(naive_results$lnZ,my_naive$lnZ)
  expect_equivalent(c(naive_results$mu),c(my_naive$mu))
  expect_equivalent(c(naive_results$alpha),c(my_naive$alpha))
  expect_equivalent(naive_results$info$iter,my_naive$iter)
  expect_equivalent(naive_results$info$max_err,my_naive$maxerr)
})



my_results <- rss_varbvsr_squarem(SiRiS = SiRiS,
                                  sigma_beta=sigb,
                                  logodds=-4.6,
                                  betahat = betahat,
                                  se = se,
                                  talpha0 = alpha_test,
                                  tmu0 = mu_test,
                                  tSiRiSr0 = SiRiSr,
                                  tolerance = 1e-4,
                                  itermax=900,
                                  verbose=T,
                                  lnz_tol=F)

test_that("SQUAREM updates are identical",{
  expect_equivalent(my_results$lnZ,mat_results$lnZ)
  expect_equivalent(my_results$mu,c(mat_results$mu))
  expect_equivalent(my_results$alpha,c(mat_results$alpha))
  expect_equivalent(my_results$iter,mat_results$info$iter)
  expect_equivalent(my_results$max_err,mat_results$info$maxerr)
})

log10oddsvec <- seq(-6,-1,0.5)
logoddsvec <- log10oddsvec*log(10)
pm <- .CallOctave('grid_rss_varbvsr_logodds',t(t(betahat)),t(t(se)),SiRiS_f,sigb,log10oddsvec,t(alpha_test),t(mu_test))

#Grid optimization has a really long parameter list
grid_options <- list(alpha=alpha_test,
                     mu=mu_test,betahat=betahat,
                     se=se,
                     SiRiS=SiRiS,
                     sigb=1,
                     logodds=logoddsvec,
                     verbose=F,
                     SiRiSr=SiRiSr,itermax=100,tolerance=1e-4,lnz_tol=F)
mr <- grid_optimize_rss_varbvsr(grid_options)
pi_mean <- marg_pi(log10odds = log10oddsvec,c(mr))
test_that("grid optimization over logodds works as expected",expect_equal(c(pi_mean),pm$pi_mean))


log10oddsvec <- seq(-3.1,-2.1,length.out = 10)
logoddsvec <- log10oddsvec*log(10)
sigb <- seq(0.8,1.2,length.out = 10)
pm <- .CallOctave('grid_rssr_varbvsr',t(t(betahat)),t(t(se)),SiRiS_f,sigb,log10oddsvec,t(alpha_test),t(mu_test))


grid_options <- list(alpha=alpha_test,
                     mu=mu_test,betahat=betahat,
                     se=se,
                     SiRiS=SiRiS,
                     sigb=sigb,
                     logodds=logoddsvec,
                     verbose=F,
                     SiRiSr=SiRiSr,itermax=100,tolerance=1e-4,lnz_tol=F)
mr_grid <- grid_optimize_rss_varbvsr(grid_options)
log_mean=rowMeans(mr_grid)
pi_mean_grid <- marg_pi(log10odds = log10oddsvec,lnz = log_mean)
sigb_mean_grid <- marg_param(colMeans(mr_grid),param = sigb)
tpi <- 10/982
test_that("2d grid optimization over sigb and logodds works as in MATLAB",expect_equal(mr_grid,pm,tolerance=1e-1))
test_that("2d grid optimization gives more or less correct values",{
  expect_equal(1,sigb_mean_grid,tolerance=0.1)
  expect_equal(10/982,c(pi_mean_grid),tolerance=0.1)})





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



