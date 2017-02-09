
rss_varbvsr_future <- function(datafile,sigb=0.058,logodds=-2.9/log(10),options=list()){
  
  stopifnot(file.exists(datafile))
  tolerance <- 1e-4
  if(is.null(options[["itermax"]])){
    itermax <- 100
  }else{
    itermax <- options[["itermax"]]
  }
  betahat_cell <- c(read_vec(datafile,"betahat"))
  se_cell <- c(read_vec(datafile,"se"))
  p <- length(betahat_cell)
  sigb_cell <- rep(sigb,p)
  logodds_cell <- rep(logodds,p)
  SiRiS_cell <- gen_SiRSi(datafile)
  if(is.null(options[["alpha"]])){
    alpha_cell <- runif(p)
    alpha_cell <- alpha_cell/sum(alpha_cell)
    mu_cell <- rnorm(p)
    options[["alpha"]] <- alpha_cell
    options[["mu"]] <- mu_cell
  }else{
    if(length(options[["alpha"]])==1){
      stopifnot(file.exists(options[["alpha"]]))
      alpha_cell <- scan(options[["alpha"]],what=numeric())
      mu_cell <- scan(options[["mu"]],what=numeric())
    }else{
      alpha_cell <- options[["alpha"]]
      mu_cell <- options[["mu"]]
    }
  }
  SiRiSr_cell <- (SiRiS_cell%*%(alpha_cell*mu_cell))@x
  return(rss_varbvsr_squarem(SiRiS = SiRiS_cell,
                             sigma_beta=sigb_cell,
                             logodds=logodds_cell,se = se_cell,
                             betahat = betahat_cell,
                             talpha0 = alpha_cell,
                             tmu0 = mu_cell,
                             tSiRiSr0 = SiRiSr_cell,
                             tolerance = tolerance))
}


rss_varbvsr_parallel_future <- function(datafiles,sigb=0.058,logodds=-2.9/log(10),options=list()){
  library(future.BatchJobs)
  library(future)
  library(h5)
  library(listenv)
  stopifnot(all(file.exists(datafiles)))
  
  if(!is.null(options[["plan"]])){
    if(options[["plan"]][["engine"]]=="PBS"){
      plan(batchjobs_torque,resources=options[["plan"]][["resources"]])
    }else{
      if(options[["plan"]][["engine"]]=="SLURM"){
        plan(batchjobs_slurm,resources=options[["plan"]][["resources"]])
      }else{
        nodes <- options[["plan"]][["resources"]][["nodes"]]
        future::plan(list(future::tweak(future::multiprocess,workers=nodes)))
      }
    }
  }else{
    future::plan(future::eager)
  }
  
  
  tolerance <- 1e-4
  if(is.null(options[["itermax"]])){
    itermax <- 100
  }else{
    itermax <- options[["itermax"]]
  }
  timevec <- numeric(length(datafiles))
  
 
  
  
  alpha_cell <- list()
  mu_cell <- list()
  
  #init_params (If we're running on the head node, we want to be as polite as possible,
  #and allow users to specify files rather than vectors)
  
  if(!is.null(options[["alpha"]])){
    if(length(options[["alpha"]])==length(datafiles)){
      alpha_cell <- options[["alpha"]]
      mu_cell <- options[["mu"]]
    }
    else{
      chrom_cell <- list()
      for(i in 1:length(datafiles)){
        chrom_cell[[i]] <- c(read_vec(datafiles[i],"chr"))
      }
      stopifnot(length(options[["alpha"]])==length(unlist(chrom_cell)),
                length(options[["mu"]])==length(unlist(chrom_cell)))
      alpha_cell <- split(options[["alpha"]],f = unlist(chrom_cell))
      mu_cell <- split(options[["mu"]],f = unlist(chrom_cell))
    }
  }
  
  resultl <- listenv()
  for(i in 1:length(datafiles)){
    cat(datafiles[i],"\n")
    if(!is.null(options[["alpha"]])){
      options[["alpha"]] <- alpha_cell[[i]]
      options[["mu"]] <- mu_cell[[i]]
      if(!is.null(options[["toFile"]])){
        cat("Submitting!\n")
        resultl[[i]] %<-% {
          cat("Batch Job Started!\n")
          tres <- rss_varbvsr_future(datafiles[i],sigb=sigb,logodds=logodds,options=options)
          outf <- h5file(options[["toFile"]][[i]],'a')
          outf["alpha"] <- tres[["alpha"]]
          outf["mu"] <- tres[["mu"]]
          outf["lnZ"] <- tres[["lnZ"]]
          outf[["iter"]] <- tres[["iter"]]
          outf[["max_err"]] <- tres[["max_err"]]
          h5close(outf)
          tres[["lnZ"]]
        }
      }else{
        resultl[[i]] %<-% rss_varbvsr_future(datafiles[i],sigb=sigb,logodds=logodds,options=options)
      }
    }else{
      resultl[[i]] %<-% rss_varbvsr_future(datafiles[i],sigb=sigb,logodds=logodds,options=options)
    }
  }
  return(future::values(resultl))
}






rss_varbvsr_parallel_squarem <- function(datafiles,sigb=0.058,logodds=-2.9/log(10),options=list()){
  
  stopifnot(all(file.exists(datafiles)))
  tolerance <- 1e-4
  if(is.null(options[["itermax"]])){
    itermax <- 100
  }else{
    itermax <- options[["itermax"]]
  }
  timevec <- numeric(length(datafiles))
  
  
  betahat_cell <- list()
  se_cell <- list()
  sigb_cell <- list()
  logodds_cell <- list()
  SiRiS_cell <- list()
  pvec <- numeric(length(datafiles))
  
  alpha_cell <- list()
  mu_cell <- list()
  SiRiSr_cell <- list()
  
  q_cell <- list()
  s_cell <- list()
  sesquare_cell <- list()  
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
    SiRiSr_cell[[i]] <- (SiRiS_cell[[i]]%*%(alpha_cell[[i]]*mu_cell[[i]]))@x
    sesquare_cell[[i]] <- se_cell[[i]]^2
    q_cell[[i]] <- betahat_cell[[i]]/sesquare_cell[[i]]
    s_cell[[i]] <- (sesquare_cell[[i]]*(sigb*sigb))/(sesquare_cell[[i]]+(sigb*sigb))
  }
  
  params_cell <- list()
  r_cell <- list()
  lnZ=-Inf
  iter <- 0
  ##Begin iteration
  
  lnZ_l <- numeric(length(datafiles))
  cat('iter   lower bound  change vars E[b] sigma2\n');
  
  for(i in 1:length(datafiles)){
    toutput <- rss_varbvsr_squarem(SiRiS = SiRiS_cell[[i]],
                                   sigma_beta=sigb_cell[[i]],
                                   logodds=logodds_cell[[i]],se = se_cell[[i]],
                                   betahat = betahat_cell[[i]],
                                   talpha0 = alpha_cell[[i]],
                                   tmu0 = mu_cell[[i]],
                                   tSiRiSr0 = SiRiSr_cell[[i]],
                                   tolerance = tolerance)
    alpha_cell[[i]] <- toutput[,1]
    mu_cell[[i]] <- toutput[,2]
    SiRiSr_cell[[i]] <- toutput[,3]
    r_cell[[i]] <-alpha_cell[[i]]*mu_cell[[i]]
    params_cell[[i]] <- c(toutput[,1],r_cell[[i]])
    
    lnZ_l[i] <- calculate_lnZ(q = q_cell[[i]],
                              r = r_cell[[i]],
                              SiRiSr = SiRiSr_cell[[i]],
                              logodds = logodds_cell[[i]],
                              sesquare = sesquare_cell[[i]],
                              alpha = alpha_cell[[i]],
                              mu = mu_cell[[i]],
                              s = s_cell[[i]],
                              sigb = sigb)
    
  }
  lnZ <- sum(lnZ_l)
  
  return(list(alpha=unlist(alpha_cell),
              mu=unlist(mu_cell),
              lnZ=lnZ,
              info=list(sigb=sigb,
                        alpha0=unlist(options[["alpha"]]),
                        mu0=unlist(options[["mu"]])),
              s=unlist(s_cell)))
}