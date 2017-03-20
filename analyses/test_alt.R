context("test alt implementation")

library(Matrix)
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
sigb_size=10
sigb_start=0.1
sigb_step=0.1
sigb <- numeric(sigb_size)
for(i in 1:sigb_size){
  sigb[i] <- sigb_start+(i-1)*sigb_step
}

logodds_size <- 10
pi_start <- 0.001
logodds_step <- 0.01
logodds <- numeric(logodds_size);
for(i in 1:logodds_size){
  pi <- pi_start+(i-1)*logodds_step
  logodds[i] <- log(pi/(1-pi))
}
tot_size <- sigb_size*logodds_size;
grid_sigb <- numeric(tot_size)
grid_logodds <- numeric(tot_size);
k <- 1
for(i in 1:sigb_size){
  for(j in 1:logodds_size){
    grid_sigb[k] <- sigb[i]
    grid_logodds[k] <- logodds[j]
    k <- k+1
  }
}

log10oddsvec <- log((10/p)/(1-(10/p)),base = 10)
logodds <- log10oddsvec*log(10)
grid_sigb <- rep(sigb,each=length(logodds))
grid_logodds <- rep(logodds,each=length(sigb))


tot_size <- length(grid_sigb)
alpha_grid <- matrix(0,nrow = p,ncol = tot_size)
mu_grid <- matrix(0,nrow = p,ncol = tot_size)
SiRiSr_grid <- matrix(0,nrow = p,ncol = tot_size)



k <- 1
for(i in 1:length(sigb)){
  for(j in 1:length(logodds)){
    t_res <- wrap_rss_varbvsr_iter(SiRiS=SiRiS,sigma_beta=sigb[i],logodds=logodds[j],
                                   betahat=betahat,
                                   se=se,
                                   alpha=alpha_test,
                                   mu=mu_test,
                                   SiRiSr=SiRiSr,
                                   reverse=F)
#    st_res <- wrap_rss_varbvsr_iter_alt(SiRiS=as.matrix(SiRiS),sigma_beta=sigb[i],logodds=logodds[j],betahat=betahat,se=se,alpha=alpha_test,mu=mu_test,SiRiSr=SiRiSr,reverse=F)
#    expect_equal(st_res,t_res)
    alpha_grid[,k] <- t_res$alpha1
    mu_grid[,k] <- t_res$mu1
    SiRiSr_grid[,k] <- t_res$SiRiSr
    k <- k+1
  }
}



t_list <- wrap_rss_varbvsr_iter_grid(SiRiS = SiRiS,sigma_beta = grid_sigb[1],logodds = grid_logodds[1],betahat = betahat,
                           se = se,alpha=alpha_test,mu = mu_test,SiRiSr = SiRiSr,reverse = F)


rss_varbvsr_naive(S)
la <- rss_varbvsr_naive(SiRiS,sigb[i],logodds[j],betahat,se,alpha_test,mu_test,SiRiSr,1e-4,100,T,F)
lb <- rss_varbvsr_alt_naive_grid(SiRiS,sigb[i],logodds[j],betahat,se,alpha_test,mu_test,SiRiSr,1e-4,100,T,F)
expect_equal(la$alpha,c(lb$alpha))
expect_equal(la$mu,c(lb$mu))
expect_equal(la$SiRiSr,c(lb$SiRiSr))
expect_equal(la$lnZ,c(lb$lnZ))
expect_equal(la$iter,lb$iter)
st_list <- list(alpha1=alpha_grid,mu1=mu_grid,SiRiSr=SiRiSr_grid)
expect_equal(alpha_grid[,1],t_list$alpha1[,1])
expect_equal(st_list,t_list)
expect_equal(t_list$alpha1,alpha_grid)
expect_equal(t_list$mu1,mu_grid)
expect_equal(t_list$SiRiSr,SiRiSr_grid)
test_that("grid search is the same between the two implementations",{
  for(i in 1:10){
#    tgrid <- rss_varbvsr_squarem(SiRiS,sigb[1],logodds[1],betahat,se,alpha_test,mu_test,SiRiSr,1e-4,100,T,F)
    grida <- grid_search_rss_varbvsr(SiRiS,sigb,logodds,betahat,se,alpha_test,mu_test,SiRiSr,1e-4,100,F,T)
  gridb <- grid_search_rss_varbvsr_alt(as.matrix(SiRiS),sigb,logodds,betahat,se,alpha_test,mu_test,SiRiSr,1e-4,100,T,F)
  gridc <- grid_search_rss_varbvsr_alt_grid(SiRiS,grid_sigb,grid_logodds,betahat,se,alpha_test,mu_test,SiRiSr,1e-4,100,T,T)
  expect_equal(grida,gridc)
  }
})

