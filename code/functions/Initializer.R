args = commandArgs(trailingOnly=TRUE)

######


chr=args[1]
dir_out=args[2]

######


########load functions#############
suppressMessages(suppressWarnings(library(rTensor)))
suppressMessages(suppressWarnings(library(MASS)))
suppressMessages(suppressWarnings(library(Matrix)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(source('./Muscle_functions.R')))
suppressMessages(suppressWarnings(library(Rcpp)))
suppressMessages(suppressWarnings(library(RSpectra)))
suppressMessages(suppressWarnings(library(qs)))
#invisible(capture.output(sourceCpp("./multiply.cpp"), type = "message"))
err<-tryCatch(invisible(capture.output(sourceCpp(paste0(dir_functions,'/multiply.cpp')), type = "message")),
         error=function(e){return(1)})
if(length(err)==1){
  
  inverse=function(x){return(chol2inv(chol(x)))}
  multiply=function(A,B){return(crossprod(t(A),B))}
  
}

R=1
print(chr)
outname=paste0(dir_out,'/data')


#For tensors
if(!(chr=="CG"|chr=="CH")){
  chr=as.numeric(args[1])
  #### Parameters ####  
  data=qs::qread(paste0(outname,chr,'_obj.qs'))
  

  
  
  btd_res=list(A=c(),B=c(),C=c())
  
  
  ntads=readRDS(paste0(dir_out,"/ntads.rds"))

  
  Lr=rep(ntads[chr],each=R)
  
  
  #####Already done###
  set.seed(1)
  hosvdres=hosvd_new(as.tensor(data),ranks=c(Lr[1],Lr[1],1))
  tmpAB=attr((ttm(hosvdres$Z,(hosvdres$U[[1]]),1)%>%ttm(.,hosvdres$U[[2]],2))[,,1],"data")
  res=eigs_sym(tmpAB, k=Lr[1])
  if(Lr[1]!=1)btd_res$A=res$vector%*%diag(res$values)
  if(Lr[1]==1)btd_res$A=as.matrix(res$vector*(res$values))
  btd_res$B=res$vector
  
  
  #%*%diag(svdres$d[1:sum(Lr)])
  saveRDS(btd_res$A,paste0(outname,chr,'_modeA_res.rds'))
  saveRDS(btd_res$B,paste0(outname,chr,'_modeB_res.rds'))
  #######################################
  
  if(Lr[1]!=1)X=as(khatri_rao_matrix2(btd_res$A,btd_res$B,R,Lr),"sparseMatrix")
  if(Lr[1]==1)X=as(khatri_rao(btd_res$A,btd_res$B),"sparseMatrix")
    #X=as(khatri_rao_matrix2(btd_res$A,btd_res$B,R,Lr),"sparseMatrix")
  Y=t(attr(rTensor::unfold(as.tensor(data),row_idx = 3,col_idx = c(2,1)) ,"data"))%>%as(.,"sparseMatrix")
  
  
  print('genertaing XTY')
  
  XTY=multiply(t(X),Y)
  
  print(paste0('saving...chr',chr,'results'))
  
  saveRDS(XTY,paste0(outname,chr,'_XTY.rds'))
  
  
  
  
  
}




#For matrix
if(chr=="CG"|chr=="CH"){
  chr=(args[1])
  #### Parameters ####  
  #CG_obj.qs
  data=qs::qread(paste0(outname,chr,'_obj.qs'))
  
 
  
  
  #####Already done###
  set.seed(1)
  svdres=svds(as(data,"sparseMatrix"),k=R)
  
  loci_loading=svdres$u*svdres$d
  saveRDS(loci_loading,paste0(outname,chr,'_Methyl_loci_res.rds'))
  #######################################
  
  
  
  
  X=as(loci_loading,"sparseMatrix")
  Y=data%>%as(.,"sparseMatrix")
  
  
  print('genertaing XTY')
  
  XTY=multiply(t(X),Y)
  
  print(paste0('saving...',chr,'results'))
  
  saveRDS(XTY,paste0(outname,chr,'_XTY.rds'))
  
  
  
  
  
}
