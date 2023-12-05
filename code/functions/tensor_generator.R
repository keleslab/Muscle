args = commandArgs(trailingOnly=TRUE)

######

chr=as.numeric(args[1])
dir_data=(args[2])
dir_functions=(args[3])
chr_num=as.numeric(args[4])
sizefile=(args[5])
binsize=as.numeric(args[6])

#ssh=args[4]
######

# dir_data=paste0(dir,"/data/example")
# dir_functions=paste0(dir,"/code/functions")
########load functions#############


setwd(dir_data)


source(paste0(dir_functions,'/Muscle_functions.R'))

pacman::p_load(abind,qs,reshape2,dplyr,rTensor,Matrix,data.table,foreach,inline,parallel,doParallel,gtools)


cl <- makePSOCKcluster(4)
registerDoParallel(cl)



chrlist=paste0("chr",c(1:(chr_num-1),"X"))


summarized_hic=qs::qread(paste0(dir_data,'/summarized_hic.qs'))

#Comeback this part


size<-fread(paste0(dir_data,'/',sizefile))
size<-size[size$V1%in%chrlist,]
size<-size[size$V1 %>% mixedorder,]


dimlist=ceiling(size$V2[size$V1%in%chrlist]/binsize)

#dimlist=readRDS('dimlist_mm9.rds')
#locpair by 3
#cell by locpair


fitted=qs::qread(paste0(dir_data,"/Matricized_HiCtensor_imputed_chr",chr,".qs"))

summarized_hic=summarized_hic[which(summarized_hic$chrom==chrlist[chr]),]



#faster tensor generator
num_cell=dim(fitted)[1]



library(foreach)
#,.export = c("sparse2dense","balance")
tmp_list=foreach(c=1:num_cell,.packages = c("Matrix",  "dplyr"))%dopar%{
  
  print(c)
  tmp=cbind(summarized_hic,fitted[c,])  
  #      tmp=hic_df %>% filter(chrom=="chr1",cell=="GM12878_1")
  tmp<-tmp%>%select(-c(chrom))
  #tmp<-tmp%>%select(-c(chrom,diag,cell))
  dimnames(tmp)[[2]]=c("i","j","val")
  tmp$i=tmp$i/binsize+1
  tmp$j=tmp$j/binsize+1
  
  aa=sparse2dense(tmp,N=dimlist[chr])
  
  return(aa) 
  
  
  
  
}



tmp_tensor=do.call("abind",c(tmp_list, list(along=3)))





print(chrlist[chr])


qs::qsave(tmp_tensor,paste0(dir_data,"/data_HiC_chr",chr,'.qs'))  

# which.chrs=paste0(c(1:chr_num),collapse =  " ")
#system(paste0('cd ',dir_functions,'; parallel -S ',ssh,'; parallel -S hodor01,hodor02,hodor03 --jobs 4 --workdir . Rscript ::: svd_logdebiased_tensor_gen_Li.R ::: ', which.chrs))
# 
# if(!is.null(ssh))system(paste0('cd ',dir_functions,'; parallel -S ',ssh,' --jobs 4 --workdir . Rscript ::: tensor_generator.R ::: ', which.chrs))
# if(is.null(ssh))system(paste0('cd ',dir_functions,'; parallel --jobs 4 --workdir . Rscript ::: tensor_generator.R ::: ', which.chrs))
# 
