
read_ccs_h5 <- function(h5filename,groupname,dataname="data",iname="ir",pname="jc"){
  require(h5)
  require(Matrix)
  h5f <- h5file(h5filename,mode = 'r')
  h5g <- h5f[groupname]
  data <- h5g[dataname][]
  i <- h5g[iname][]
  p <- h5g[pname][]
  h5close(h5f)
  return(sparseMatrix(i=i+1,p=p,x=data))
}


read_vec <- function(h5filename,datapath){
  require(h5)
  h5f <- h5file(h5filename,'r')
  data <- h5f[datapath][]
  h5close(h5f)
  return(data)
}

gen_SiRSi <- function(h5filename){
  require(h5)
  require(Matrix)
  spmat <- read_ccs_h5(h5filename,"R")
  se <- read_vec(h5filename,"se")
  return(SiRSi(spmat,1/se))
}