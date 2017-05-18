
#' Run Regression with Summary Statistics over a grid of hyperparamter values
#' 
#' This is an optimized implementation of RSS. Given a correlation matrix
#'  and univariate summary statistics, it will return a `data.frame` with the 
#'  variational lower bound evaluated over a grid of hyperparameter values.
#'  
#'  @param R Correlation matrix (`p` by `p`)
#'    This can be a standard dense R matrix, or a sparse matrix from the package `Matrix`
#' @param sigma_beta Hyperparameter sa is the prior variance of regression 
#' coefficients for variables that are included in the model. If it is not provided as input, then it is set
#' to an evenly spaced grid (log-scale) of values from 10^0.
#' 
#' @param logodds Hyperparameter logodds is the prior log-odds that a variable is included in the regression model;
#'  it is defined as logodds = log(q/(1-q)), where q is the prior probability that a variable is included in the
#'  regression model. Note that we use the natural logarithm instead of the base-10 logarithm.
#'  If logodds is not provided as input, then it is set to an evenly spaced grid of values from `-14` to `-4.6`
#' @param betahat a length p vector specifying the univariate regression coefficients from the association study.
#' @param se a length p vector specifying the univariate standard errors from the association study


rssr_grid <- function(R,betahat,se,sigma_beta=NULL,logodds=NULL,options=list(use_squarem=T,
                                                                        use_lnztol=T,
                                                                        tolerance=1e-3,
                                                                        itermax=200,
                                                                        parallel=F,
                                                                        nthreads=1)){
  p <- lenghth(betahat)
  stopifnot(length(se)==p,
            nrow(R)==ncol(R),
            nrow(R)==p,
            all(se>0),tolerance>0,
            itermax>=0)
  
  
  if(typeof(R)=="S4"){
    isSparse <- T
  }else{
    isSparse<- F
  }
  if(isSparse){
    siris <- SiRSi(R,1/se)
  }else{
    siris <- SiRSi_d(R,1/se)
  }
  
  if(is.null(sigma_beta)){
    sigma_beta <- 10^seq(-2,0.5,length.out = 10)
  }
  if(is.null(logodds)){
    pi <- 10^seq(-6,-2,length.out=10)
    logodds <- log(pi/(1-pi))
  }
  if(length(logodds)!=length(sigma_beta)){
    paramdf <- list(logodds=logodds,sigma_beta=sigma_beta) %>% cross_d() %>% dplyr::distinct()
    logodds <- paramdf$logodds
    sigma_beta <- paramdf$sigma_beta
  }
  
  if(is.null(options[["alpha0"]])){
    alpha0 <- ralpha(p)
  }else{
    alpha0 <- options[["alpha0"]]
  }
  if(is.null(options[["mu0"]])){
    mu0 <- rmu(p)
  }else{
    mu0 <- options[["mu0"]]
  }
  sirisr=c(siris%*%(alpha0*mu0))
  if(is.null(options[["tolerance"]])){
    tolerance <- 1e-3
  }else{
    tolerance <- options[["tolerance"]]
  }
  if(is.null(options[["use_squarem"]])){
    use_squarem <- T
  }else{
    use_squarem <- options[["use_squarem"]]
  }
  if(is.null(options[["use_lnztol"]])){
    use_lnztol<- T
  }else{
    use_lnztol <- options[["use_lnztol"]]
  }
  if(is.null(options[["itermax"]])){
    itermax <- 200
  }else{
    itermax <- options[["itermax"]] 
  }
  
  
  if(!isSparse){
    if(use_squarem){
      result <- grid_search_rss_varbvsr(SiRiS = siris,sigma_beta = sigma_beta,logodds = logodds,betahat = betahat,
                                        se = ,talpha0 = alpha0,tmu0 = mu0,tSiRiSr0 = sirisr,
                                        tolerance = tolerance,itermax = itermax,verbose = F,lnz_tol = use_lnztol)
    }else{
      result  <-grid_search_rss_varbvsr_naive(SiRiS = siris,sigma_beta = sigma_beta,logodds = logodds,betahat = betahat,
                                              se = ,talpha0 = alpha0,tmu0 = mu0,tSiRiSr0 = sirisr,
                                              tolerance = tolerance,itermax = itermax,verbose = F,lnz_tol = use_lnztol)
    }
  }else{
    if(use_squarem){
      result <- grid_search_rss_varbvsr_sp(SiRiS = siris,sigma_beta = sigma_beta,logodds = logodds,betahat = betahat,
                                        se = ,talpha0 = alpha0,tmu0 = mu0,tSiRiSr0 = sirisr,
                                        tolerance = tolerance,itermax = itermax,verbose = F,lnz_tol = use_lnztol)
    }else{
      result  <-grid_search_rss_varbvsr_naive_sp(SiRiS = siris,sigma_beta = sigma_beta,logodds = logodds,betahat = betahat,
                                              se = ,talpha0 = alpha0,tmu0 = mu0,tSiRiSr0 = sirisr,
                                              tolerance = tolerance,itermax = itermax,verbose = F,lnz_tol = use_lnztol)
    }
  }
  return(result)
}

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




 


