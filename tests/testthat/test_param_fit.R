context("Fitting parameters")

data(betahat)

data(se)

data(alpha_test)
data(mu_test)
data(R_shrink)

betahat <- c(betahat)
se <- c(se)
mu_test <- c(mu_test)
alpha_test <- c(alpha_test)
t_SiRiS <- SiRSi(R_shrink,1/se)
SiRiS_f <- as.matrix(SiRSi(R_shrink,1/se))
p <- length(betahat)
SiRiSr=c(SiRiS_f%*%(alpha_test*mu_test))
SiRiS  <-as(SiRiS_f,"dgCMatrix")


test_that("We can estimate a known value of pi by grid search",{
  sigb <- 1
  log10oddsvec <- seq(-6,-1,0.5)
  logodds <- log10oddsvec*log(10)
  # input_h5files <- dir("/media/nwknoblauch/Data/GTEx/1kg_LD",full.names = T)[1]
  # output_h5files <- gsub("1kg_LD","1kg_IBD",input_h5files)
  # file.remove(output_h5files)
  
  
  options <- list(plan=list(engine="MC",resources=list(nodes=2)),
                  betahat=betahat,
                  se=se,
                  logodds=logodds,
                  alpha=alpha_test,
                  mu=mu_test,
                  SiRiS=SiRiS,
                  sigb=sigb,
                  datafile=tempfile(),
                  tolerance=1e-4)
  datafiles <- tempfile()
  mmat <- rss_varbvsr_parallel_future(options=options)
  lnzmat <- gen_lnzmat(mmat[[1]],logodds ,sigb)
  
  
  expect_equal(c(marg_pi(log10oddsvec,lnzmat)),(10/p),tolerance=0.01)
}
)

test_that("We can estimate a value of sigb by grid search",{
  sigb <- seq(0.5,1.5,0.1)
  log10oddsvec <- log((10/p)/(1-(10/p)),base = 10)
  logodds <- log10oddsvec*log(10)
  options <- list(plan=list(engine="MC",resources=list(nodes=2)),
                  betahat=betahat,
                  se=se,
                  logodds=logodds,
                  alpha=alpha_test,
                  mu=mu_test,
                  SiRiS=SiRiS,
                  sigb=sigb,
                  datafile=tempfile(),
                  tolerance=1e-4)
  datafiles <- tempfile()
  mmat <- rss_varbvsr_parallel_future(options=options)
  lnzmat <- gen_lnzmat(mmat[[1]],logodds ,sigb)
  gsigb <- marg_param(c(lnzmat),param = sigb)

  expect_equal(gsigb,1,tolerance=0.1)
})





test_that("parallel and serial grid search work the same in 2d",{
  
  
  sigb <- seq(0.5,1.5,0.1)
  log10oddsvec <- seq(-6,-1,0.5)
  logodds <- log10oddsvec*log(10)
  serial_mat <- matrix(0,nrow = length(logodds),ncol = length(sigb))
  for(i in 1:nrow(serial_mat)){
    for(j in 1:ncol(serial_mat)){
      serial_mat[i,j] <- rss_varbvsr_squarem_iter(SiRiS,sigb[j],logodds[i],betahat,se,alpha_test,mu_test,SiRiSr,1e-3,100,T)
    }
  }
  
  p_df <- grid_search_rss_varbvsr(SiRiS = SiRiS,sigma_beta = sigb,logodds = logodds,
                                  betahat = betahat,se = se,talpha0 = alpha_test,tmu0 = mu_test,
                                  tSiRiSr0 = SiRiSr,tolerance = 1e-3,itermax = 100,verbose = F,lnz_tol = T)
  
  
  
  expect_equal(c(serial_mat),c(p_df$lnZ))
  expect_equal(rep(logodds,length(sigb)),p_df$logodds)
  expect_equal(sort(rep(sigb,length(logodds))),p_df$sigb)
  
})



test_that("parallel and serial grid search work the same in 1d (sigb)",{
  
  sigb <- seq(0.5,1.5,0.1)
  log10oddsvec <- log((10/p)/(1-(10/p)),base = 10)
  logodds <- log10oddsvec*log(10)
  serial_mat <- matrix(0,nrow = length(logodds),ncol = length(sigb))
  for(i in 1:nrow(serial_mat)){
    for(j in 1:ncol(serial_mat)){
      serial_mat[i,j] <- rss_varbvsr_squarem_iter(SiRiS,sigb[j],logodds[i],betahat,se,alpha_test,mu_test,SiRiSr,1e-3,100,T)
    }
  }
  
  p_df <- grid_search_rss_varbvsr(SiRiS = SiRiS,sigma_beta = sigb,logodds = logodds,
                                  betahat = betahat,se = se,talpha0 = alpha_test,tmu0 = mu_test,
                                  tSiRiSr0 = SiRiSr,tolerance = 1e-3,itermax = 100,verbose = F,lnz_tol = T)
  
  
  
  expect_equal(c(serial_mat),c(p_df$lnZ))
  expect_equal(rep(sigb,length(logodds)),p_df$sigb)
  expect_equal(rep(logodds,length(sigb)),p_df$logodds)
})


test_that("parallel and serial grid search work the same in 1d (logodds)",{
  
  
  sigb <- c(1)
  log10oddsvec <- seq(-6,-1,0.5)
  logodds <- log10oddsvec*log(10)
  serial_mat <- matrix(0,nrow = length(logodds),ncol = length(sigb))
  for(i in 1:nrow(serial_mat)){
    for(j in 1:ncol(serial_mat)){
      serial_mat[i,j] <- rss_varbvsr_squarem_iter(SiRiS,sigb[j],logodds[i],betahat,se,alpha_test,mu_test,SiRiSr,1e-3,100,T)
    }
  }
  p_df <- grid_search_rss_varbvsr(SiRiS = SiRiS,sigma_beta = sigb,logodds = logodds,
                                  betahat = betahat,se = se,talpha0 = alpha_test,tmu0 = mu_test,
                                  tSiRiSr0 = SiRiSr,tolerance = 1e-3,itermax = 100,verbose = F,lnz_tol = T)
  
  
  
  expect_equal(c(serial_mat),p_df$lnZ)
  expect_equal(rep(sigb,length(logodds)),p_df$sigb)
  expect_equal(rep(logodds,length(sigb)),p_df$logodds)
})



