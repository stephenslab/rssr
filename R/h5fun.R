

  
#' read a sparse matrix from HDF5
#' read a sparse matrix from HDF5 stored in Compressed Column Storage (CCS) format
#' @template h5fun
#' @export 
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
  require(h5)
  if(check){
    h5f <- h5file(h5filename,'r')
    if(existsGroup(h5f,"SiRiS")){
      cat("Found!\n")
      h5close(h5f)
      return(read_ccs_h5(h5filename,"SiRiS",isSymmetric=FALSE))
    }
    h5close(h5f)
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