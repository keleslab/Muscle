
source("config_file_model.R")

##########################################################


pacman::p_load(RColorBrewer,MASS,Matrix,dplyr,cluster,rTensor,reshape2,Rcpp,RcppArmadillo,foreach,inline,parallel,doParallel,RSpectra,qs)
pacman::p_load(Rcpp,RSpectra,qs,RColorBrewer)





system(paste0("mkdir -p ",dir_out))

#Rank of HiC tensors
Rank=as.numeric(readRDS(paste0(dir_data,'/Rank.rds')))
# number of chromosomes(human:23, mouse:20)
chr_num=as.numeric(readRDS(paste0(dir_data,'/chr_num.rds')))



source(paste0(dir_functions,'/Muscle_functions.R'))

#invisible(capture.output(sourceCpp(paste0(dir_functions,'/multiply.cpp')), type = "message"))

#err<-tryCatch(invisible(capture.output(sourceCpp(paste0(dir_functions,'/multiply.cpp')), type = "message")),
#         error=function(e){return(1)})
#if(length(err)==1){
  
#  inverse=function(x){return(chol2inv(chol(x)))}
#  multiply=function(A,B){return(crossprod(t(A),B))}
  
#}


if(Bulk_exist==TRUE){
  
#  LR=c()
#  for (i in c(1:22, "X")) {
    
    #bulk_tad=
#    LR=c(LR,dim(bulk_tad)[1])}
#  saveRDS(LR,paste0(dir_out,'/ntads.rds'))
LR=readRDS(paste0(dir_out,'/ntads.rds'))  
  
}


#normalize HiC and Methyl


print("Reading in data...")

if(modality=="HiC"){
    datatensor=mapply(function(x){
    tmp=log(qs::qread(paste0(dir_data,"/data_HiC_chr",x,'.qs'))+0.0001)
    return(tmp)},as.list(1:chr_num))
    if(Bulk_exist==FALSE){
      dimlist=lapply(datatensor, function(x)dim(x)[1]) %>% unlist
      #LR=floor(dimlist/Rank)
       LR=floor(dimlist/10)
      saveRDS(LR,paste0(dir_out,'/ntads.rds'))
      
    }

}


if(modality=="All"|modality=="HiC+CG"){
  
  datatensor=mapply(function(x){
    tmp=qs::qread(paste0(dir_data,"/data_HiC_chr",x,'.qs'));
    
    nonzeroind=Reduce("union",apply(tmp,3,function(x)which(rowMeans(x!=0)!=0)))
    tmp=tmp[nonzeroind,nonzeroind,]
    
    tmp=log(tmp+0.0001)
    return(array(apply(tmp,3,scalemat),dim=dim(tmp)))},as.list(1:chr_num))
  if(Bulk_exist==FALSE){
    dimlist=lapply(datatensor, function(x)dim(x)[1]) %>% unlist
       LR=floor(dimlist/10)
    saveRDS(LR,paste0(dir_out,'/ntads.rds'))
    
  }
  

  
  
  tmp=qs::qread(paste0(dir_data,"/data_methyl_","CG",'.qs'))
  tmp=naomit_matrix(tmp)
  datatensor[['CG']]=apply(tmp,2,scalemat)
  
if(Rank_mCG==0){  
  svd_res=RSpectra::svds(datatensor[['CG']],k = exploration_rank)
  nn_sv=svd_res$d[svd_res$d>0]
  (nn_sv ) %>% log %>% plot
  
  pdf(file = paste0(dir_data,"/svd_mCG_plot.pdf"),width = 7,height = 7) 
  (nn_sv ) %>% log %>% plot(.,main="singular value plot of mCG",ylab="log(lambda)",xlab="Rank")
  abline(v=10,col=2)
  abline(v=20,col=3)
  abline(v=30,col=4)
  abline(v=40,col=5)
  legend('topright', legend=c("Rank=10", "Rank=20", "Rank=30", "Rank=40"),
         col=c(2,3,4,5),lty = 1,cex=1)
  dev.off()
  
  
  

  cat("Please specify the mCG matrix rank value based on the singular value 'svd_mCG_plot.pdf' and hit enter (skipping will give rank=30) : ")
  Rank_mCG <- as.numeric(readLines(con="stdin", 1))
  cat(Rank_mCG, "\n")        
  if(is.na(Rank_mCG)){Rank_mCG=30}                   
}
  


  
  
  
  
if(modality=="All"){
  tmp=qs::qread(paste0(dir_data,"/data_methyl_","CH",'.qs'))
  tmp=naomit_matrix(tmp)
  datatensor[['CH']]=apply(tmp,2,scalemat)
  
if(Rank_mCH==0){  
  svd_res=RSpectra::svds(datatensor[['CH']],k = exploration_rank)
  nn_sv=svd_res$d[svd_res$d>0]
  (nn_sv ) %>% log %>% plot
  
  pdf(file = paste0(dir_data,"/svd_mCH_plot.pdf"),width = 7,height = 7) 
  (nn_sv ) %>% log %>% plot(.,main="singular value plot of mCH",ylab="log(lambda)",xlab="Rank")
  abline(v=10,col=2)
  abline(v=20,col=3)
  abline(v=30,col=4)
  abline(v=40,col=5)
  legend('topright', legend=c("Rank=10", "Rank=20", "Rank=30", "Rank=40"),
         col=c(2,3,4,5),lty = 1,cex=1)
  dev.off()
  
  
  
  
  cat("Please specify the mCH matrix rank value based on the singular value 'svd_mCH_plot.pdf' and hit enter (skipping will give rank=30) : ")
  Rank_mCH  <- as.numeric(readLines(con="stdin", 1))
  cat(Rank_mCH, "\n")

  if(is.na(Rank_mCH)){Rank_mCH=30}
   } 
    
    
    
  Rank=ceiling(mean(c(Rank,Rank_mCG,Rank_mCH)))
  
  
}
  
  if(modality=="HiC+CG"){Rank=ceiling(mean(c(Rank,Rank_mCG)))}  
  
}



Muscle=function(data,R,tol=0.001,maxiter=5,dir_out,dir_functions,chr_num=20,modality="HiC",ssh=NULL){
  
  
  outname=paste0(dir_out,'/data')
  
  
  
  for(k in 1:R){
    
    if(k==1){
      
      if(modality=="All"){obj=data;mapply(function(x,y){qs::qsave(x,paste0(outname,y,'_obj.qs'))},obj,as.list(c(1:chr_num,'CG','CH')))}
      if(modality=="HiC+CG"){obj=data;mapply(function(x,y){qs::qsave(x,paste0(outname,y,'_obj.qs'))},obj,as.list(c(1:chr_num,'CG')))}
      if(modality=="HiC"){obj=data;mapply(function(x,y){qs::qsave(x,paste0(outname,y,'_obj.qs'))},obj,as.list(c(1:chr_num)))}
  
      rm(data)
      }
    print(paste0("Module ",k," start"))
    
    #rank 1 update
    tmp=rankone_Muscle(data=obj,k,LR=LR,tol=tol,maxiter=maxiter,dir_out=dir_out,dir_functions=dir_functions,chr_num,modality,ssh,GNU)
    
    if(k==1){btd_res=tmp$btd_res;
    if(modality=="All"|modality=="HiC+CG"){
      methyl_res_CG=tmp$methyl_res_CG
      if(modality=="All"){methyl_res_CH=tmp$methyl_res_CH}
      
      } 
             }
    
    if(k!=1){for(chr in 1:chr_num){
      btd_res[[chr]]$A=cbind(btd_res[[chr]]$A,tmp$btd_res[[chr]]$A)
      btd_res[[chr]]$B=cbind(btd_res[[chr]]$B,tmp$btd_res[[chr]]$B)
      btd_res[[chr]]$C=cbind(btd_res[[chr]]$C,tmp$btd_res[[chr]]$C)
    };
      if(modality=="All"|modality=="HiC+CG"){  
      methyl_res_CG$Methyl_loci_loading=cbind(methyl_res_CG$Methyl_loci_loading,tmp$methyl_res_CG$Methyl_loci_loading)
      methyl_res_CG$Methyl_cell_loading=cbind(methyl_res_CG$Methyl_cell_loading,tmp$methyl_res_CG$Methyl_cell_loading)
      if(modality=="All"){methyl_res_CH$Methyl_loci_loading=cbind(methyl_res_CH$Methyl_loci_loading,tmp$methyl_res_CH$Methyl_loci_loading)
        methyl_res_CH$Methyl_cell_loading=cbind(methyl_res_CH$Methyl_cell_loading,tmp$methyl_res_CH$Methyl_cell_loading)}
      }
    }
    saveRDS(btd_res,paste0(dir_out,'/HiC_result_final_','_rank',k,'.rds'))
    if(modality=="All"|modality=="HiC+CG"){  
    saveRDS(methyl_res_CG,paste0(dir_out,'/methyl_CG_result_final_','_rank',k,'.rds'))
    if(modality=="All"){saveRDS(methyl_res_CH,paste0(dir_out,'/methyl_CH_result_final_','_rank',k,'.rds'))}
    }
    
    
    if(modality=="All"){obj=mapply(function(x,y,z){tmp=x-y;qs::qsave(tmp,paste0(outname,z,'_obj.qs'));return(tmp)},obj,tmp$Theta,as.list(c(1:chr_num,'CG','CH')))}
    if(modality=="HiC+CG"){obj=mapply(function(x,y,z){tmp=x-y;qs::qsave(tmp,paste0(outname,z,'_obj.qs'));return(tmp)},obj,tmp$Theta,as.list(c(1:chr_num,'CG')))}
    if(modality=="HiC"){obj=mapply(function(x,y,z){tmp=x-y;qs::qsave(tmp,paste0(outname,z,'_obj.qs'));return(tmp)},obj,tmp$Theta,as.list(c(1:chr_num)))}
    rm(tmp)
    
    
  }
  
  print("Muscle over saving results...")
  saveRDS(btd_res,paste0(dir_out,'/HiC_result_final_','_rank',k,'.rds'))
  if(modality=="All"|modality=="HiC+CG"){  
  saveRDS(methyl_res_CG,paste0(dir_out,'/methyl_CG_result_final_','_rank',k,'.rds'))
  if(modality=="All"){saveRDS(methyl_res_CH,paste0(dir_out,'/methyl_CH_result_final_','_rank',k,'.rds'))    }
  
  }
  
  
  if(modality=="HiC")print(paste0("Results can be found at ",dir_out,'/HiC_result_final_','_rank',Rank,'.rds'))
  if(modality=="HiC+CG")print(paste0("Results can be found at ",dir_out,'/HiC_result_final_','_rank',Rank,'.rds',"  and  ",dir_out,'/methyl_CG_result_final_','_rank',Rank,'.rds'))
  if(modality=="All")print(paste0("Results can be found at ",dir_out,'/HiC_result_final_','_rank',Rank,'.rds',"  and  ",dir_out,'/methyl_CG_result_final_','_rank',Rank,'.rds'," and ",dir_out,'/methyl_CH_result_final_','_rank',Rank,'.rds'))

} 





set.seed(1)

Muscle(datatensor,R=Rank,tol=tol,maxiter=maxiter,dir_out=dir_out,dir_functions = dir_functions,chr_num=chr_num,modality = modality,ssh)
