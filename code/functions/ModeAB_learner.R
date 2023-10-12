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

if(!(chr=="CG"|chr=="CH")){
  chr=as.numeric(args[1])
  data=qs::qread(paste0(outname,chr,'_obj.qs'))
  
  btd_res=readRDS(paste0(outname,chr,'_btd_res.rds'))
  
  
  ntads=readRDS(paste0(dir_out,"/ntads.rds"))
  Lr=rep(ntads[chr],each=R)
  
  #Lr=rep(round(dim(data)[1]/10),R)
  
  cY=((as.tensor(data) %>% ttm(.,t(as.matrix(btd_res$C)),3) %>% attr(.,"data")) [,,1]) %>% as.matrix
  
  res=eigs_sym(cY, k=Lr[1])
  #btd_res$A=res$vector%*%diag(res$values)
  if(Lr[1]!=1)btd_res$A=res$vector%*%diag(res$values)
  if(Lr[1]==1)btd_res$A=as.matrix(res$vector*(res$values))
  btd_res$B=as.matrix(res$vector)
  
  print(paste0('saving...chr',chr,'results'))
  saveRDS(btd_res$A,paste0(outname,chr,'_modeA_res.rds'))
  saveRDS(btd_res$B,paste0(outname,chr,'_modeB_res.rds'))
  
  ##########################################
  
  
  if(Lr[1]!=1)X=as(khatri_rao_matrix2(btd_res$A,btd_res$B,R,Lr),"sparseMatrix")
  if(Lr[1]==1)X=as(khatri_rao(btd_res$A,btd_res$B),"sparseMatrix")
  Y=t(attr(rTensor::unfold(as.tensor(data),row_idx = 3,col_idx = c(2,1)) ,"data"))%>%as(.,"sparseMatrix")
  
  
  print('genertaing XTY')
  
  XTY=multiply(t(X),Y)
  
  print(paste0('saving...chr',chr,'results'))
  
  saveRDS(XTY,paste0(outname,chr,'_XTY.rds'))
  saveRDS(as.numeric(multiply(t(X),X)),paste0(outname,chr,'_XTX.rds'))
  
  
  
  
  
  
}

if(chr=="CG"|chr=="CH"){
  chr=(args[1])
  data=qs::qread(paste0(outname,chr,'_obj.qs'))
  
  
  methyl_res=readRDS(paste0(outname,'CG','_Methyl_res.rds'))
  
  
  Yc=(multiply(as(data,"sparseMatrix"),as(methyl_res$Methyl_cell_loading,"sparseMatrix"))) %>% as.matrix
  

  loci_loading=Yc
  saveRDS(loci_loading,paste0(outname,chr,'_Methyl_loci_res.rds'))
  #######################################################
  
  
  X=as(loci_loading,"sparseMatrix")
  Y=data%>%as(.,"sparseMatrix")
  
  
  print('genertaing XTY')
  
  XTY=multiply(t(X),Y)
  
  print(paste0('saving...',chr,'results'))
  
  saveRDS(XTY,paste0(outname,chr,'_XTY.rds'))
  
  
  
  
  
}
