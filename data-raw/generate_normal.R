library(RSSReQTL)
library(devtools)



rssr_all_norm <- function(R,betahat_mat,se_mat,sigb,tolerance=1e-3,lnz_tol=T,itermax=200,fgeneid=NULL,n=1,use_squarem=T){
  
  library(dplyr)
  
  library(progress)
  library(BBmisc)
  library(purrr)
  
  # paramdf <- list(sigb=sigb,logodds=logodds) %>% cross_d() %>% distinct()
  # sigb <- paramdf$sigb
  # logodds <- paramdf$logodds  
  
  stopifnot(ncol(R)==nrow(betahat_mat),
            ncol(betahat_mat)==ncol(se_mat))
  
  if(is.null(fgeneid)){
    fgeneid <- 1:ncol(betahat_mat)
  }
  
  
  num_sigb <- length(sigb)
  ng <- ncol(betahat_mat)
  
  retl <- list()
  p <- nrow(betahat_mat)
  
  
  pb <- progress_bar$new(total=ng)
  retdfl <- list()
  nind=n
  for(i in 1:ng){
    SiRiS <- SiRSi_d(R,Si=1/se_mat[,i])
    # SiRiS <- SiRSi_d(R,se = se_mat[,i],betahat = betahat_mat[,i],n = nind)
    
    alpha0 <- ralpha(p = p)
    mu0 <-rmu(p)
    SiRiSr <- (SiRiS%*%(alpha0*mu0))
    mfgeneid <- fgeneid[i]
    if(use_squarem){
      retdfl[[i]] <- grid_search_rss_varbvsr_norm(SiRiS = SiRiS,
                                                  sigma_beta = sigb,
                                                  betahat = betahat_mat[,i],
                                                  se = se_mat[,i],talpha0 = alpha0,
                                                  tmu0 = mu0,tSiRiSr0 = SiRiSr,
                                                  tolerance = tolerance,itermax=itermax,verbose = F,
                                                  lnz_tol = lnz_tol) %>% 
        mutate(pve=pve/nind,fgeneid=mfgeneid)
    }else{
      retdfl[[i]] <- grid_search_rss_varbvsr_naive_norm(SiRiS = SiRiS,
                                                        sigma_beta = sigb,
                                                        betahat = betahat_mat[,i],
                                                        se = se_mat[,i],talpha0 = alpha0,
                                                        tmu0 = mu0,tSiRiSr0 = SiRiSr,
                                                        tolerance = tolerance,itermax=itermax,verbose = F,
                                                        lnz_tol = lnz_tol) %>% 
        mutate(pve=pve/nind,fgeneid=mfgeneid)
      
      
    }
    pb$tick()
  }
  return(bind_rows(retdfl))
  #    data_frame(lnZ=lnZvec,pi=pivec,alpha_mean=alpha_meanvec,sigb=sigbvec,fgeneid=fgeneidvec,pve=pvevec))
}





simulate_estimate <- function(betamat,tparam_df,R,useNorm=T,n){
  p <- ncol(R)
  SNP <- MASS::mvrnorm(n=n,Sigma = R,mu = rep(0,p),empirical = F)
  residvec <- calc_residvec(tparam_df,SNP,betamat)
  
  residmat <- sim_residmat(n=n,residvec = residvec)
  
  ymat <- scale(SNP%*%betamat+residmat,center=T,scale = F)
  
  betahat_mat <- map_beta_exp(SNP,ymat)
  
  se_mat <- map_se_exp(SNP,ymat,betahat_mat)
  
  bounds_sigb <- c(min(tparam_df$tsigb),max(tparam_df$tsigb))
  bounds_logodds <- c(min(tparam_df$tlogodds),max(tparam_df$tlogodds))
  
  if(useNorm){
    paramdf <- list(sigb=seq(bounds_sigb[1]/2,bounds_sigb[2]*1.5,length.out = 100))
    rss_res <- rssr_all_norm(R=R,
                             betahat_mat = betahat_mat,
                             se_mat = se_mat,
                             sigb = paramdf$sigb,
                             tolerance = 1e-3,
                             lnz_tol = T,
                             itermax = 3000,n = n,use_squarem = T)
  }else{
    paramdf <- list(sigb=seq(bounds_sigb[1]/2,bounds_sigb[2]*1.5,length.out = 10),
                    logodds=seq(bounds_logodds[1]-0.3,bounds_logodds[2]+0.3,length.out=10)) %>% cross_d() %>% distinct()
    rss_res <- rssr_all(R=R,
                        betahat_mat = betahat_mat,
                        se_mat = se_mat,
                        sigb = paramdf$sigb,logodds = paramdf$logodds,
                        tolerance = 1e-3,
                        lnz_tol = T,
                        itermax = 3000,
                        n = n,
                        use_squarem = T)
  }
  
  rss_summ <-group_by(rss_res,fgeneid) %>% mutate(w=normalizeLogWeights(lnZ)) %>% 
    summarise(mean_alpha=sum(w*alpha_mean),
              mean_pve=sum(w*pve),
              mean_sigb=sum(w*sigb)) %>% ungroup() %>% mutate(isNorm=useNorm,n=n) %>% inner_join(tparam_df,by="fgeneid")
  return(rss_summ)
  
  
}


library(purrr)
library(dplyr)
data("R_shrink")
n <- 500
p <- nrow(R_shrink)

pve.seq_norm <- c(.8)
sigb.seq_norm <- as.numeric(seq(0.3,0.9,length.out = 10))
pi.seq_norm <- 1


nreps <- 1
fparams_norm <- list(tpve=pve.seq_norm,tsigb=sigb.seq_norm,tpi=pi.seq_norm) %>% cross_d() %>% distinct()

tparam_df_norm <- bind_rows(replicate(nreps,fparams_norm,simplify = F)) %>% 
  group_by(tpve,tsigb,tpi) %>% 
  mutate(replicate=1:n()) %>% 
  ungroup() %>% 
  mutate(fgeneid=1:n(),tlogodds=Inf)


betamat_norm <- sim_betamat(tparam_df_norm,p)
n <- 382  
R_sd <- as.matrix(R_shrink)
SNP <- scale(MASS::mvrnorm(n=n,mu=rep(0,p),Sigma=R_shrink),center=T,scale=F)


residvec_SNP_norm <- calc_residvec(tparam_df_norm,SNP,betamat_norm)

residmat_SNP_norm <- sim_residmat(n=n,residvec = residvec_SNP_norm)

ymat_SNP_norm <- scale(SNP%*%betamat_norm+residmat_SNP_norm,center=T,scale = F)


betahat_mat_norm <- map_beta_exp(SNP,ymat_SNP_norm)

se_mat_norm <- map_se_exp(SNP,ymat_SNP_norm,betahat_mat_norm)


devtools::use_data(betahat_mat_norm)
devtools::use_data(se_mat_norm)
devtools::use_data(tparam_df_norm)





