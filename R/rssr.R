


# rssr <- function(R,betahat_mat,se_mat)

rss_varbvsr <- function(options=list()){
  stopifnot(length(options[["sigb"]])==1,
            length(options[["logodds"]])==1,
            !is.null(options[["SiRiS"]]),
            !is.null(options[["betahat"]]),
            !is.null(options[["alpha"]]),
            !is.null(options[["mu"]]),
            !is.null(options[["se"]]),
            length(options[["se"]])==length(options[["betahat"]]),
            length(options[["betahat"]])==length(options[["mu"]]))
  if(options[["method"]]=="naive"){  
    run_time <- system.time(int_res <- rss_varbvsr_naive(SiRiS = options[["SiRiS"]],
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
    
  }else{
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
            !is.null(options[["SiRiSr"]]),
            !is.null(options[["itermax"]]),
            !is.null(options[["lnz_tol"]]),
            length(options[["se"]])==length(options[["betahat"]]),
            length(options[["betahat"]])==length(options[["mu"]]))
  
  sigbb <- options[["sigb"]]
  logoddsb <- options[["logodds"]]
  fsigb <- runif(1,min = sigbb[1],max = sigbb[2])
  flogodds <- runif(1,min =logoddsb[1],max = logoddsb[2])
  parv=c(flogodds,fsigb)
  moptim <- optim(par = parv,fn = wrap_rss_varbvs_squarem_optim,
                  lower=c(logoddsb[1],sigbb[1]),upper=c(logoddsb[2],sigbb[2]),
                  SiRiS=options[["SiRiS"]],
                  betahat = options[["betahat"]],
                  se = options[["se"]],
                  talpha0 = options[["alpha"]],
                  tmu0 = options[["mu"]],
                  tSiRiSr0 = options[["SiRiSr"]],
                  tolerance = options[["tolerance"]],
                  itermax=options[["itermax"]],
                  lnz_tol = options[["lnz_tol"]],method="L-BFGS-B")
  return(moptim)
}




 


