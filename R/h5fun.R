

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

  
#' read a sparse matrix from HDF5
#' read a sparse matrix from HDF5 stored in Compressed Column Storage (CCS) format
#' @template h5fun

read_ccs_h5 <- function(h5filename,groupname,dataname="data",iname="ir",pname="jc",isSymmetric=NULL){
  require("h5")
  require("Matrix")
  h5f <- h5::h5file(h5filename,mode = 'r')
  stopifnot(existsGroup(h5f,groupname)|existsDataSet(h5f,groupname))
  h5g <- h5f[groupname]
  h5_attrs <- list.attributes(h5g)
  if(!is.null(dim(h5g))){
    spmat <- as(Matrix(data = h5g[],sparse = T,forceCheck = F),"dgCMatrix")
    h5close(h5g)
    return(spmat)
  }
  if("isSymmetric" %in% h5_attrs){
    storeSymmetric <- h5attr(h5g,"isSymmetric")=="TRUE"    
  }else{
    storeSymmetric <- F
  }
  if("Layout" %in% h5_attrs){
    layout <- h5attr(h5g,"Layout")
    stopifnot(layout=="CCS")
  }
  data <- h5g[dataname][]
  i <- h5g[iname][]
  p <- h5g[pname][]
  h5::h5close(h5f)
  if("Dimensions" %in% h5_attrs){
    dims <- h5attr(h5g,"Dimensions")
    if(storeSymmetric){
      if(!is.null(isSymmetric)){
        if(!isSymmetric){
          spmat <- Matrix::sparseMatrix(i=i,p=p,x=data,index1 = F,dims = dims)
          return(genSymm(spmat))
        }
      }
      return(Matrix::sparseMatrix(i=i,p=p,x=data,index1 = F,dims = dims,symmetric = T))
    }
    return(Matrix::sparseMatrix(i=i,p=p,x=data,index1 = F,dims=dims,symmetric = F))
  }else{
    return(Matrix::sparseMatrix(i=i,p=p,x=data,index1 = F,symmetric = F))
  }
}

read_lnzmat_h5 <- function(h5filename){
  require(h5)
  stopifnot(file.exists(h5filename))
  h5f <- h5file(h5filename,'r')
  
  
}




#' write a sparse matrix to HDF5
#' write a sparse matrix to HDF5 stored in the Compressed Column Storage (CCS) format)
#' @template h5fun
write_ccs_h5 <- function(h5filename,spmat,groupname,dataname="data",iname="ir",pname="jc",compression_level=8,symmetric=F){
  require(h5)
  h5f <- h5::h5file(h5filename,'a')
  h5g <- createGroup(h5f,groupname)
  #Matrix::sparseMatrix(i=i+1,p=p,x=data))
  i <- spmat@i
  p <- spmat@p
  data <- spmat@x
  h5::createAttribute(.Object = h5g,attributename = "Layout",data="CCS")
  h5::createAttribute(.Object = h5g,attributename = "Dimensions",data=dim(spmat))
  h5::createAttribute(.Object = h5g,attributename = "isSymmetric",data=ifelse(symmetric,"TRUE","FALSE"))
  id <- h5::createDataSet(.Object = h5g,
                          datasetname = iname,
                          data=i,
                          chunksize=as.integer(length(i)/10),
                          maxdimensions = NA_integer_,
                          compression=as.integer(compression_level))
  pd <- h5::createDataSet(.Object = h5g,
                          datasetname = pname,
                          data=p,
                          chunksize=as.integer(length(p)/10),
                          maxdimensions = NA_integer_,
                          compression=as.integer(compression_level))
  
  dd <- h5::createDataSet(.Object = h5g,
                          datasetname = dataname,
                          data=data,
                          chunksize=as.integer(length(data)/10),
                          maxdimensions = NA_integer_,
                          compression=as.integer(compression_level))
  
  h5close(h5f)
  return(T)
}

#' @describeIn read_ccs_h5 Read vector from an HDF5 file
read_vec <- function(h5filename,datapath){
  requireNamespace("h5")
  h5f <- h5::h5file(h5filename,'r')
  data <- h5f[datapath][]
  h5::h5close(h5f)
  return(data)
}

#' @describeIn read_ccs_h5 Read vector from an HDF5 file
write_vec <- function(h5filename,vec,datapath){
  requireNamespace("h5")
  h5f <- h5::h5file(h5filename,'a')
  data <- h5f[datapath] <- vec
  h5::h5close(h5f)
}



#' @describeIn read_ccs_h5 Read in R and standard error, and generate the SiRiS matrix
gen_SiRSi <- function(h5filename,Rpath=NULL,se=NULL,check=T){
  requireNamespace("h5")
  if(check){
    h5f <- h5::h5file(h5filename,'r')
    if(h5::existsGroup(h5f,"SiRiS")){
      cat("Found!\n")
      h5::h5close(h5f)
      return(read_ccs_h5(h5filename,"SiRiS",isSymmetric=FALSE))
    }
    h5::h5close(h5f)
  }
  if(is.null(Rpath)){
    Rpath <- "R"
  }
  spmat <- read_ccs_h5(h5filename,Rpath,isSymmetric = FALSE)
  if(is.null(se)){
    se <- read_vec(h5filename,"se")
  }
  return(SiRSi(spmat,1/se))
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



