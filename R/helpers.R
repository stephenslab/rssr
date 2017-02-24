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


marg_pi_mat <- function(lnzmat,log10odds=NULL,sigb=NULL){
  
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
    

#'Flexible function for specifying data and hyperparameters for RSS
#'@param datafile path to HDF5 formatted file with betahat, se, and R
#'@param options list containing the remaining parameters that have already been specified.
#'@details The elements of options are :
#'betahat: vector of length p pre-specifying betahat (in this case, don't read betahat from datafile)
#'se: vector of length p pre-specifying se (in this case, don't read se from datafile)
#'R: A sparse matrix of class dgCMatrix, specifying LD (note, this is distinct from the SiRiS matrix)
#'SiRiS: A sparse matrix of class dgCMatrix, specifying the SiRiS matrix
#'tolerance: the tolerance for convergence of the VB algorithm
#'maxiter: the maximum number of iterations to perform
#'alphafile file containing initial values of alpha (txt file)
#'mu file containing initial values of alpha (txt file)
#'alpha: a length p vector specifying initial values of alpha
#'mu: a length p vector specifying initial values of mu
#'Rfile: HDF5 file to use instead of datafile for LD matrix
#'betahatfile: HDF5 file to use instead of datafile for betahat vector
#'sefile: HDF5 file to use instead of datafile for se vector
#'Rpath: path in HDF5 file to look for R matrix (may be name of group in case of sparse matrix or name of data in case of dense matrix)
#'betahatpath: path in HDF5 file to look for betahat vector 
#'sepath: path in HDF5 file to look for se vector 

prep_rss <- function(datafile=NULL,options=list(),chunk=NULL,tot_chunks=NULL){
  require(h5)
  options[["datafile"]] <- prep_list(options[["datafile"]],datafile,chunk,tot_chunks)
  options[["tolerance"]] <- prep_list(options[["tolerance"]],1e-4,chunk,tot_chunks)
  options[["itermax"]] <- prep_list(options[["itermax"]],100,chunk,tot_chunks)
  
  options[["betahatfile"]] <- prep_list(options[["betahatfile"]],options[["datafile"]],chunk,tot_chunks)
  options[["betahatpath"]] <- prep_list(options[["betahatpath"]],"betahat",chunk,tot_chunks)
  options[["betahat"]] <- prep_list(options[["betahat"]],  c(read_vec(options[["betahatfile"]],options[["betahatpath"]])),chunk,tot_chunks)
  
  options[["sefile"]] <- prep_list(options[["sefile"]],options[["datafile"]],chunk,tot_chunks)
  options[["sepath"]] <- prep_list(options[["sepath"]],"se",chunk,tot_chunks)
  options[["se"]] <- prep_list(options[["se"]],  c(read_vec(options[["sefile"]],options[["sepath"]])),chunk,tot_chunks)
  options[["verbose"]] <- prep_list(options[["verbose"]],T,chunk,tot_chunks)
  options[["lnz_tol"]] <- prep_list(options[["lnz_tol"]],T,chunk,tot_chunks)
  p <- length(options[["se"]])
  stopifnot(length(options[["se"]])==length(options[["betahat"]]))
  if(is.null(options[["SiRiS"]])){
    options[["Rfile"]] <- prep_list(options[["Rfile"]],options[["datafile"]],chunk,tot_chunks)
    options[["Rpath"]] <- prep_list(options[["Rpath"]],"R",chunk,tot_chunks)
    options[["se"]] <- prep_list(options[["se"]],  c(read_vec(options[["sefile"]],options[["sepath"]])),chunk,tot_chunks)
    options[["SiRiS"]] <- gen_SiRSi(options[["Rfile"]],options[["Rpath"]],options[["se"]],check=F)
  }
  options[["alphafile"]] <- prep_list(options[["alphafile"]],NULL,chunk,tot_chunks)
  options[["mufile"]] <- prep_list(options[["mufile"]],NULL,chunk,tot_chunks)
  
  if(!is.null(options[["alphafile"]])){
    stopifnot(file.exists(options[["alphafile"]]))
    options[["alpha"]] <- scan(options[["alphafile"]],what=numeric())
  }else{
    options[["alpha"]] <- prep_list(options[["alpha"]],NULL,chunk,tot_chunks)
    if(is.null(options[["alpha"]])){
      options[["alpha"]] <- ralpha(p)
    }
  }
  if(!is.null(options[["mufile"]])){
    stopifnot(file.exists(options[["mufile"]]))
    options[["mu"]] <- scan(options[["mufile"]],what=numeric())
  }else{
    options[["mu"]] <- prep_list(options[["mu"]],NULL,chunk,tot_chunks)
    if(is.null(options[["mu"]])){
      options[["mu"]] <- rmu(p)
    }
  }
  if(is.null(options[["SiRiSr"]])){
    options[["SiRiSr"]] <- (options[["SiRiS"]]%*%(options[["alpha"]]*options[["mu"]]))@x
  }
  if(is.null(options[["sigb"]])){
    options[["sigb"]] <- 0.058
  }
  if(is.null(options[["logodds"]])){
    options[["logodds"]] <- -2.9/log(10)
  }
  
  # tempf <- tempfile()
  # th <- h5file(tempf,'a')
  # thd <- createDataSet(th,datasetname = "tdn",data=1:5,chunksize=2L,compression=4L,maxdimensions = NA_integer_)
  # h5close(th)
  # file.remove(tempf)
  
  
  
  return(options)
}

frac_resolved <- function(job_list){
  require(future)
  totsize <- 0
  tot_resolved <- 0
  for(i in 1:length(job_list)){
    for(j in 1:length(job_list[[i]])){
      for(k in 1:length(job_list[[i]][[j]])){
        if(resolved(job_list[[i]][[j]][[k]])){
          tot_resolved <- tot_resolved+1
        }
        totsize <- totsize+1
      }
    }
  }
  return(tot_resolved/totsize)
}



