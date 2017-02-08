
#' read a sparse matrix from HDF5
#' read a sparse matrix from HDF5 stored in Compressed Column Storage (CCS) format
#' @template h5fun
#' @export 
read_ccs_h5 <- function(h5filename,groupname,dataname="data",iname="ir",pname="jc"){
  requireNamespace("h5")
  requireNamespace("Matrix")
  h5f <- h5::h5file(h5filename,mode = 'r')
  h5g <- h5f[groupname]
  data <- h5g[dataname][]
  i <- h5g[iname][]
  p <- h5g[pname][]
  h5::h5close(h5f)
  return(Matrix::sparseMatrix(i=i+1,p=p,x=data))
}

#' write a sparse matrix to HDF5
#' write a sparse matrix to HDF5 stored in the Compressed Column Storage (CCS) format)
#' @template h5fun

write_ccs_h5 <- function(h5filename,spmat,groupname,dataname="data",iname="ir",pname="jc",compression_level=8){
  require(h5)
  h5f <- h5::h5file(h5filename,'a')
  h5g <- createGroup(h5f,groupname)
  #Matrix::sparseMatrix(i=i+1,p=p,x=data))
  i <- spmat@i
  p <- spmat@p
  data <- spmat@x
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
  return(data)
}

#' @describeIn read_ccs_h5 Read in R and standard error, and generate the SiRiS matrix
gen_SiRSi <- function(h5filename,check=T){
  requireNamespace("h5")
  require(h5)
  if(check){
    h5f <- h5file(h5filename,'r')
    if(existsGroup(h5f,"SiRiSi")){
      h5close(h5f)
      return(read_ccs_h5(h5filename,"SiRiSi"))
    }
    h5close(h5f)
  }
  spmat <- read_ccs_h5(h5filename,"R")
  se <- read_vec(h5filename,"se")
  return(SiRSi(spmat,1/se))
}