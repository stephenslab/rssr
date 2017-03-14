#' Convenience function to call two RSS updates
#' @template rssr
#' @param o_mat a p by 3 matrix consisting of the concatenation of the three vectors: alpha0,mu0 and SiRiSr0, the initial values of alpha,mu and SiRiSr.
squarem_update <- function(SiRiS,
                           sigma_beta,
                           logodds,
                           betahat,
                           se,o_mat,reverse){
  alpha0 <- o_mat[,1]
  mu0 <- o_mat[,2]
  SiRiSr0 <- o_mat[,3]
  resmat <- rss_varbvsr_iter_naive_reference(SiRiS  =SiRiS,
                        sigma_beta = sigma_beta,
                        logodds = logodds,
                        betahat = betahat,
                        se = se,
                        alpha = alpha0,
                        mu = mu0,
                        SiRiSr = SiRiSr0,reverse)
  
  alpha1 <- resmat[,1]
  mu1 <- resmat[,2]
  SiRiSr1 <- resmat[,3]
  
  sresmat <- rss_varbvsr_iter_naive_reference(SiRiS=SiRiS,
                         sigma_beta=sigma_beta,
                         logodds = logodds,
                         betahat=betahat,
                         se = se,
                         alpha = alpha1,
                         mu = mu1,
                         SiRiSr = SiRiSr1,
                         reverse)
  return(array(c(resmat,sresmat),c(dim(resmat),2)))
}





#' RSS
#' Main function for running RSS in serial (with SQUAREM updates)
#' @param datafiles a vector of HDF5 files with one file per chromosome.(See details for constraints on the layout of these files)
#' @param sigb The prior SD of the regression coefficients (if included), scalar
#' @param logodds  logodds: the prior log-odds (i.e. log(prior PIP/(1-prior PIP))) of inclusion for each SNP, p by 1
#' @param options a list containing optional arguments (including initial values for mu and alpha)
#' @return A list with the following elements:
#' lnZ: scalar, the variational lower bound of the marginal log likelihood (up to some constant)
#'       alpha: p by 1, variational estimates of the posterior inclusion probabilities 
#'      mu: p by 1, posterior means of the additive effects (given snp included)
#'      s: p by 1, posterior variances of the additive effects (given snp included)
#'       info: list with following fields 
#'               - iter: integer, number of iterations
#'               - maxerr: the maximum relative difference between the parameters at the last two iterations
#'               - sigb: scalar, the maximum likelihood estimate of sigma_beta
#'               - loglik: iter by 1, the variational lower bound at each iteration
rss_varbvsr_bigmem_squarem <- function(datafiles,sigb=0.058,logodds=-2.9/log(10),options=list()){
  
  stopifnot(all(file.exists(datafiles)))
  tolerance <- 1e-4
  if(is.null(options[["itermax"]])){
    itermax <- 200
  }else{
    itermax <- options[["itermax"]]
  }
  timevec <- numeric(length(datafiles))
  
  alpha_r_norm2 <- numeric(length(datafiles))
  alpha_v_norm2 <- numeric(length(datafiles))
  mu_r_norm2 <- numeric(length(datafiles))  
  mu_v_norm2 <- numeric(length(datafiles))
  ntimevec <- numeric(length(datafiles))
  timevec3 <- numeric(length(datafiles))
  maxerr_uni <- numeric(length(datafiles))
  absr_uni <- numeric(length(datafiles))
  asum_uni<- numeric(length(datafiles))
  loglik <- numeric(0)
  
  
  
  betahat_cell <- list()
  se_cell <- list()
  sigb_cell <- list()
  logodds_cell <- list()
  SiRiS_cell <- list()
  pvec <- numeric(length(datafiles))
  
  alpha_cell <- list()
  mu_cell <- list()
  SiRiSr_cell <- list()
  params_cell <- list()
  
  q_cell <- list()
  s_cell <- list()
  sesquare_cell <- list()  
  chrom_cell <- list()
  for(i in 1:length(datafiles)){
    # cat(datafiles[i],"\n")
    betahat_cell[[i]] <- c(read_vec(datafiles[i],"betahat"))
    se_cell[[i]] <- c(read_vec(datafiles[i],"se"))
    pvec[i] <- length(betahat_cell[[i]])
    sigb_cell[[i]] <- rep(sigb,pvec[i])
    logodds_cell[[i]] <- rep(logodds,pvec[i])
    SiRiS_cell[[i]] <- gen_SiRSi(datafiles[i])
    chrom_cell[[i]] <- c(read_vec(datafiles[i],"chr"))
  }
  
  
  #init_params
  if(is.null(options[["alpha"]])){
    for(i in 1:length(datafiles)){
      alpha_cell[[i]] <- runif(pvec[i])
      alpha_cell[[i]] <- alpha_cell[[i]]/sum(alpha_cell[[i]])
      mu_cell[[i]] <- rnorm(pvec[i])
    }
    options[["alpha"]] <- unlist(alpha_cell)
    options[["mu"]] <- unlist(mu_cell)
  }else{
    alpha_cell <- split(options[["alpha"]],f = unlist(chrom_cell))
    mu_cell <- split(options[["mu"]],f = unlist(chrom_cell))
  }
  for(i in 1:length(datafiles)){
    params_cell[[i]] <- c(alpha_cell[[i]],alpha_cell[[i]]*mu_cell[[i]])
    SiRiSr_cell[[i]] <- (SiRiS_cell[[i]]%*%(alpha_cell[[i]]*mu_cell[[i]]))@x
    sesquare_cell[[i]] <- se_cell[[i]]^2
    q_cell[[i]] <- betahat_cell[[i]]/sesquare_cell[[i]]
    s_cell[[i]] <- (sesquare_cell[[i]]*(sigb*sigb))/(sesquare_cell[[i]]+(sigb*sigb))
  }
  
  
  lnZ=-Inf
  alpha_r <- list()
  alpha_v <- list()
  mu_r <- list()
  mu_v <- list()
  
  r_list <- list()
  
  alpha_tmp <- list()
  mu_tmp <- list()
  SiRiSr_tmp <- list()
  lnZ_l <- numeric(length(datafiles))
  
  alpha0_cell <- list()
  mu0_cell <- list()
  SiRiSr0_cell <- list()
  params0_cell <- list()
  iter <- 0
  ##Begin iteration
  
  
  cat('iter   lower bound  change vars E[b] sigma2\n');
  while(iter<itermax){
    
    lnZ0 <- lnZ
    
    for(i in 1:length(datafiles)){
      alpha0_cell[[i]] <- alpha_cell[[i]]
      mu0_cell[[i]] <- mu_cell[[i]]
      params0_cell[[i]] <- params_cell[[i]]
      SiRiSr0_cell[[i]] <- SiRiSr_cell[[i]]
    }
    reverse <- iter%%2!=0    
    
    for(i in 1:length(datafiles)){
      # cat(i,"\n")
      gc()
      
      ntimevec[i] <- system.time(toutput <- squarem_update(SiRiS=SiRiS_cell[[i]],
                                                           sigma_beta=sigb_cell[[i]],
                                                           logodds=logodds_cell[[i]],se = se_cell[[i]],
                                                           betahat = betahat_cell[[i]],
                                                           o_mat = cbind(alpha_cell[[i]],mu_cell[[i]],SiRiSr_cell[[i]]),reverse))["elapsed"]
      alpha_0 <- alpha_cell[[i]]
      alpha_1 <- toutput[,1,1]
      alpha_2 <- toutput[,1,2]
      
      mu_0 <- mu_cell[[i]]
      mu_1 <- toutput[,2,1]
      mu_2 <- toutput[,2,2]
      
      alpha_r[[i]] <- alpha_1-alpha_0
      mu_r[[i]] <- mu_1-mu_0
      
      alpha_v[[i]] <- (alpha_2-alpha_1)-alpha_r[[i]]
      mu_v[[i]] <- (mu_2-mu_1)-mu_r[[i]]
      
      alpha_r_norm2[i] <- sum(alpha_r[[i]]^2)
      mu_r_norm2[i] <- sum(mu_r[[i]]^2)
      
      alpha_v_norm2[i] <- sum(alpha_v[[i]]^2)
      mu_v_norm2[i] <- sum(mu_v[[i]]^2)
      
      alpha_tmp[[i]] <- alpha_2
      mu_tmp[[i]] <- mu_2
      SiRiSr_tmp[[i]] <- toutput[,3,2]
    }
    
    
    mtp <- -sqrt(sum(alpha_r_norm2)+sum(mu_r_norm2))/sqrt(sum(alpha_v_norm2)+sum(mu_v_norm2))
    
    if(mtp >=-1 ){
      for(i in 1:length(datafiles)){
        2+2
      }
    } else{
      for(i in 1:length(datafiles)){
        alpha_tmp[[i]] <- alpha_cell[[i]]-2*mtp*alpha_r[[i]]+(mtp^2)*alpha_v[[i]]
        mu_tmp[[i]] <- mu_cell[[i]]-2*mtp*mu_r[[i]]+(mtp^2)*mu_v[[i]]
        SiRiSr_tmp[[i]] <- (SiRiS_cell[[i]] %*% (alpha_tmp[[i]]*mu_tmp[[i]]))@x
      }
    }
    
    for(i in 1:length(datafiles)){
      # cat(i,"\n")
      timevec3[i]<-system.time(toutput<-rss_varbvsr_iter_naive_reference(SiRiS=SiRiS_cell[[i]],
                                                    sigma_beta=sigb_cell[[i]],
                                                    logodds=logodds_cell[[i]],
                                                    betahat = betahat_cell[[i]],se=se_cell[[i]],
                                                    alpha0 = alpha_tmp[[i]],mu0 = mu_tmp[[i]],
                                                    SiRiSr0 = SiRiSr_tmp[[i]],reverse))["elapsed"]
      
      alpha_cell[[i]] <- toutput[,1]
      mu_cell[[i]] <- toutput[,2]
      SiRiSr_cell[[i]] <- toutput[,3]  
      
      
      r_list[[i]] <- toutput[,1]*toutput[,2]
      params_cell[[i]] <- c(toutput[,1],r_list[[i]])
      lnZ_l[i] <- calculate_lnZ(q = q_cell[[i]],
                                r = r_list[[i]],
                                SiRiSr = SiRiSr_cell[[i]],
                                logodds = logodds_cell[[i]],
                                sesquare = sesquare_cell[[i]],
                                alpha = alpha_cell[[i]],
                                mu = mu_cell[[i]],
                                s = s_cell[[i]],
                                sigb = sigb)
      J_tmp <- which(params_cell[[i]]>1e-6)
      err_tmp <- relerr(params_cell[[i]][J_tmp],params0_cell[[i]][J_tmp])
      maxerr_uni[i] <- max(err_tmp)
    }
    
    
    lnZ <- sum(lnZ_l)
    
    if((mtp<(-1))&(lnZ<lnZ0)){
      num_bt=0
      while((lnZ<lnZ0)&(num_bt<10)){
        mtp <- 0.5*(mtp-1)
        for(i in 1:length(datafiles)){
          
          alpha_tmp3 <- alpha0_cell[[i]] - 2*mtp*alpha_r[[i]]+(mtp^2)*alpha_v[[i]]
          mu_tmp3 <- mu0_cell[[i]]-2*mtp*mu_r[[i]]+(mtp^2)*mu_v[[i]]
          SiRiSr_tmp3 <- (SiRiS_cell[[i]]%*%(alpha_tmp3*mu_tmp3))@x
          
          loutput <- rss_varbvsr_iter_naive_reference(SiRiS=SiRiS_cell[[i]],
                                 sigma_beta=sigb_cell[[i]],
                                 logodds=logodds_cell[[i]],
                                 betahat = betahat_cell[[i]],se=se_cell[[i]],
                                 alpha0 = alpha_tmp3,mu0 = mu_tmp3,SiRiSr0 = SiRiSr_tmp3,reverse)
          alpha_cell[[i]] <- loutput[,1]
          mu_cell[[i]] <- loutput[,2]
          r <- loutput[,1]*loutput[,2]
          SiRiSr_cell[[i]] <- loutput[,3]
          lnZ_l[i] <- calculate_lnZ(q = q_cell[[i]],
                                    r = r_list[[i]],
                                    SiRiSr = SiRiSr_cell[[i]],
                                    logodds = logodds_cell[[i]],
                                    sesquare = sesquare_cell[[i]],
                                    alpha = alpha_cell[[i]],
                                    mu = mu_cell[[i]],
                                    s = s_cell[[i]],
                                    sigb = sigb)
          
        }
        lnZ <- sum(lnZ_l)
        num_bt <- num_bt+1
      }
    }
    for(i in 1:length(datafiles)){
      r=alpha_cell[[i]]*mu_cell[[i]]
      params_cell[[i]] <- c(alpha_cell[[i]],r)
      J_tmp <- which(params_cell[[i]]>1e-6)
      err_tmp <- relerr(params_cell[[i]][J_tmp],params0_cell[[i]][J_tmp])
      maxerr_uni[i] <- max(err_tmp)
      absr_uni[i] <- max(abs(r))
      asum_uni[i] <- sum(alpha_cell[[i]])
    }
    
    maxerr <- max(maxerr_uni)
    absr <- max(absr_uni) 
    asum <- round(sum(asum_uni))
    
    loglik <- c(loglik,lnZ)
    if(lnZ<lnZ0){
      cat('\n');
      cat('WARNING: the log variational lower bound decreased by %+0.2e\n',lnZ0-lnZ);
      final_mu <- unlist(mu0_cell)
      lnZ <- lnZ0
      break
    }else {
      if(maxerr<tolerance){
        final_alpha <- unlist(alpha_cell)
        final_mu <- unlist(mu_cell)
        cat('\n');
        cat('Convergence reached: maximum relative error %+0.2e\n',maxerr);
        cat('The log variational lower bound of the last step increased by %+0.2e\n',lnZ-lnZ0);
        break
      }
    }
    iter=iter+1
    if(iter==itermax){
      final_alpha <- unlist(alpha_cell)
      final_mu <- unlist(mu_cell)
      break
    }
    status = sprintf('%4d %+13.6e %0.1e %4d %0.2f %5.2f\n',iter,lnZ,maxerr,asum,absr,sigb*sigb);
    cat(status);
  }
  return(list(alpha=unlist(alpha_cell),
              mu=unlist(mu_cell),
              lnZ=lnZ,
              info=list(iter=iter,
                        loglik=loglik,
                        sigb=sigb,
                        alpha0=unlist(options[["alpha"]]),
                        mu0=unlist(options[["mu"]])
              ),
              s=unlist(s_cell)))
}


