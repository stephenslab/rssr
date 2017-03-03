
rss_varbvsr_future <- function(options=list()){
  if(options[["verbose"]]){
    cat('iter   lower bound  change vars E[b] sigma2\n');
  }
  stopifnot(length(options[["sigb"]])==1,
            length(options[["logodds"]])==1,
            !is.null(options[["SiRiS"]]),
            !is.null(options[["betahat"]]),
            !is.null(options[["alpha"]]),
            !is.null(options[["mu"]]),
            !is.null(options[["se"]]),
            length(options[["se"]])==length(options[["betahat"]]),
            length(options[["betahat"]])==length(options[["mu"]]))
            
  run_time <- system.time(int_res <- rss_varbvsr_squarem(SiRiS = options[["SiRiS"]],
                                                         sigma_beta=options[["sigb"]],
                                                         logodds=options[["logodds"]],
                                                         betahat = options[["betahat"]],
                                                         se = options[["se"]],
                                                         talpha0 = options[["alpha"]],
                                                         tmu0 = options[["mu"]],
                                                         tSiRiSr0 = options[["SiRiSr"]],
                                                         tolerance = options[["tolerance"]],
                                                         itermax=options[["itermax"]],
                                                         verbose=options[["verbose"]],
                                                         lnz_tol = options[["lnz_tol"]]))
  int_res[["run_time"]] <- run_time
  return(int_res)
}

rss_varbvsr_optim <- function(options=list()){
  library(stats)
  stopifnot(length(options[["sigb"]])==2,
            length(options[["logodds"]])==2,
            !is.null(options[["SiRiS"]]),
            !is.null(options[["betahat"]]),
            !is.null(options[["alpha"]]),
            !is.null(options[["mu"]]),
            !is.null(options[["se"]]),
            length(options[["se"]])==length(options[["betahat"]]),
            length(options[["betahat"]])==length(options[["mu"]]))
  
  
  

}



rss_varbvsr_parallel_grid <- function(options=list()){
  
  resultl <- list()
  datafile <- options[["datafile"]]
  for(i in 1:length(options[["datafile"]])){
    data_opts <- prep_rss(datafile[i],options=options,chunk=i,tot_chunks=length(datafile))
    resultl[[i]] <- grid_search_rss_varbvsr(SiRiS = data_opts[["SiRiS"]],
                                sigma_beta=unlist(data_opts[["sigb"]]),
                                logodds=unlist(data_opts[["logodds"]]),
                                betahat = data_opts[["betahat"]],
                                se = data_opts[["se"]],
                                talpha0 = data_opts[["alpha"]],
                                tmu0 = data_opts[["mu"]],
                                tSiRiSr0 = data_opts[["SiRiSr"]],
                                tolerance = data_opts[["tolerance"]],
                                itermax=data_opts[["itermax"]],
                                verbose=data_opts[["verbose"]],
                                lnz_tol = data_opts[["lnz_tol"]])
  }
  return(resultl)
  
  
}

grid_optimize_rss_varbvsr <- function(options=list()){
  lnzmat <- matrix(data=NA,nrow=length(options[["logodds"]]),ncol = length(options[["sigb"]]))
  for(i in 1:nrow(lnzmat)){
    for(j in 1:ncol(lnzmat)){
      cat("logodds:",options[["logodds"]][i],"\n")
      cat("sigb:",options[["sigb"]][j],"\n")
      int_res <- rss_varbvsr_squarem(SiRiS = options[["SiRiS"]],
                                     sigma_beta=options[["sigb"]][j],
                                     logodds=options[["logodds"]][i],
                                     betahat = options[["betahat"]],
                                     se = options[["se"]],
                                     talpha0 = options[["alpha"]],
                                     tmu0 = options[["mu"]],
                                     tSiRiSr0 = options[["SiRiSr"]],
                                     tolerance = options[["tolerance"]],
                                     itermax=options[["itermax"]],
                                     verbose=options[["verbose"]],
                                     lnz_tol = options[["lnz_tol"]])
      lnzmat[i,j] <- int_res[["lnZ"]]
    }
  }
  return(lnzmat)
}
  

 

rss_varbvsr_parallel_future <- function(options=list()){
  
  library(future.BatchJobs)
  library(future)
  library(h5)
  
  if(!is.null(options[["plan"]])){
    if(options[["plan"]][["engine"]]=="PBS"){
      cat("PBS/Torque detected\n")
      plan(batchjobs_torque,resources=options[["plan"]][["resources"]])
    }else{
      if(options[["plan"]][["engine"]]=="SLURM"){
        cat("SLURM detected\n")
        plan(batchjobs_slurm,resources=options[["plan"]][["resources"]])
      }else{
        if(options[["plan"]][["engine"]]=="MC"){
          cat("Multisession parallelism\n")
          nodes <- as.integer(options[["plan"]][["resources"]][["nodes"]])
          future::plan(multiprocess, workers = nodes)
        }else{
          cat("Custon parallelism :",options[["plan"]][["engine"]],"\n")
          tplan <- options[["plan"]][["plan"]]
          future::plan(tplan)
        }
        #        future::plan(list(future::tweak(future::multiprocess,workers=nodes)))
      }
    }
  }else{
    cat("Defaulting to serial execution\n")
    future::plan(future::eager)
  }
  
  
  #init_params (If we're running on the head node, we want to be as polite as possible,
  #and allow users to specify files rather than vectors)
  #  chunk_opts <- prep_rss_chunks(datafiles = datafiles,logoddsvec = logodds,sigbvec = sigb,options = options)
  stopifnot((!is.null(options[["logodds"]])),
            (!is.null(options[["sigb"]])))
  resultl <- list()
  datafiles <- unlist(options[["datafile"]])
  sigbvec <- options[["sigb"]]
  logoddsvec <- options[["logodds"]]
  for(i in 1:length(datafiles)){
    cat("File: ",i,"of ",length(datafiles),"\n")
    resultl[[i]] <- list()
    for(j in 1:length(logoddsvec)){
      resultl[[i]][[j]] <- list()
      options[["logodds"]] <- logoddsvec[j]
      for(k in 1:length(sigbvec)){
        options[["sigb"]]<-sigbvec[k]
        cat("logodds: ",j,"of ",length(logoddsvec),"\n")
        cat("sigb: ",k,"of ",length(sigbvec),"\n")
        resultl[[i]][[j]][[k]] <- future({
          data_opts <- prep_rss(datafiles[i],options=options,chunk=i,tot_chunks=length(datafiles))
          tres <- rss_varbvsr_future(options = data_opts)
          tres[["lnZ"]]
        })
      }
    }
  }
  
  cat("Waiting on Results")
  fr <- frac_resolved(resultl)
  cat(fr,"\n")
  while(fr<1){
    cat("Waiting on:",(1-fr)*100,"% of  results \n")
    Sys.sleep(5)
    fr <- frac_resolved(resultl)
  }
  cat("All Resolved!\n")
  return(resultl)
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