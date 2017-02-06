relerr <- function(x1, x2){
  if (is.null(x1) | is.null(x2)){
    return(0)  
  }else{
    return(abs(x1-x2)/(abs(x1)+abs(x2)+.Machine$double.eps))
  }
}




read_option_vec <- function(optionvec,pvec){
  npvec <- cumsum(c(1,pvec[-length(pvec)]))
  cpvec <- cumsum(pvec)
  optionlist <- list()
  for(i in 1:length(pvec)){
    optionlist[[i]] <- optionvec[npvec[i]:cpvec[i]]
    stopifnot(length(optionlist[[i]])==pvec[i])
  }
  return(optionvec)
}


#' Function to create an HDF5 dataset
#' @param h5filename output file name
#' @template rssr
#' 
create_h5_data <- function(h5filename,R,betahat,se,chr=NULL){
  cat("Writing file ",h5filename,"\n")
  if(is.null(chr)){
    chr <- sample(1:10000,1)
  }
  stopifnot(nrow(R)==length(betahat),length(betahat)==length(se))
  #make sure that chromosome vector is same length as other data vectors
  chr <- chr[rep(1,length(betahat))]
  write_vec(h5filename,chr,"chr")
  write_vec(h5filename,betahat,"betahat")
  write_vec(h5filename,se,"se")
  write_ccs_h5(h5filename=h5filename,spmat = R,groupname = "R")
  cat("Finished writing file ",h5filename,"\n")
  return(T)
}