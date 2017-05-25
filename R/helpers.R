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
    

libpath <- function() {
  cat(sprintf(
    "%s/rssr/libs/rssr%s",
    installed.packages()["rssr","LibPath"][1],
    .Platform$dynlib.ext
  ))
}
