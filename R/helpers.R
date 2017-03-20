relerr <- function(x1, x2){
  if (is.null(x1) | is.null(x2)){
    return(0)  
  }else{
    return(abs(x1-x2)/(abs(x1)+abs(x2)+.Machine$double.eps))
  }
}


#' Normalize log weights
#' @param lnZ unnormalized log probability
#' normalizelogweights takes as input an array of unnormalized
# log-probabilities logw and returns normalized probabilities such
# that the sum is equal to 1.
normalizeLogWeights <- function (lnZ,na.rm=T) {
  
  # Guard against underflow or overflow by adjusting the
  # log-probabilities so that the largest probability is 1.
  const <- max(lnZ,na.rm = na.rm)
  w <- exp(lnZ - const)
  
  # Normalize the probabilities.
  return(w/sum(w,na.rm = na.rm))
}

marg_pi<- function(log10odds=NULL,lnz=NULL){
  log10odds <- log10odds[!is.na(lnz)]
  lnz <- lnz[!is.na(lnz)]
  lw=normalizeLogWeights(lnz)
  pi=1/(1+10^(-log10odds))
  return(lw%*%pi)
}



marg_param <- function(lnZ,param){
  param <- param[!is.na(lnZ)]
  lnZ <- lnZ[!is.na(lnZ)]
  normw <- normalizeLogWeights(lnZ)
  mean_param <- c(normw%*%param)
  return(c(mean_param))
}
 


gen_lnzmat <- function(matlist,logoddsvec,sigbvec){
  require(future)
  lnzmat <- matrix(0,nrow = length(matlist),ncol = length(matlist[[1]]))
  for(i in 1:length(logoddsvec)){
    for(j in 1:length(sigbvec)){
      tres <- value(matlist[[i]][[j]])
      if(!is.null(tres)){
        lnzmat[i,j] <- tres
      }else{
        lnzmat[i,j] <- NA
      }
    }
  }
  return(lnzmat)
}

  


ralpha <- function(p){
  alpha <- runif(p)
  return(alpha/sum(alpha))
}

rmu <- function(p){
  return(rnorm(p))
}
  

#'Helper function  to specify data and hyperparameters for the scenario 
#where LD,betahat,and se are in the same files but split into chunks,
#where alpha and mu, if available, are in separate txt files
prep_list <- function(x,default,chunk=NULL,chunk_max=NULL){
  if(!is.null(x)){
    if(!is.null(chunk)){
      if(length(x)==chunk_max){
        if(typeof(x)=="list")
          return(x[[chunk]])
      }
    }
    return(x)
  }
  return(default)
}
    




# gen_test_data <- function(){
#   data("betahat")
#   data("se")
# }

libpath <- function() {
  cat(sprintf(
    "%s/rssr/libs/rssr%s",
    installed.packages()["rssr","LibPath"][1],
    .Platform$dynlib.ext
  ))
}
