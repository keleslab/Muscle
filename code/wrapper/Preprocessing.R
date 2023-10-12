

source("config_file_preprocess.R")


#######

setwd(dir_data)

saveRDS(chr_num,paste0(dir_data,'/chr_num.rds'))

pacman::p_load(RSpectra,qs,reshape2,dplyr,rTensor,Matrix,data.table,gtools)
options(scipen = 4)


##################################


cat("Reading in scHi-C data... \n")

hic_df=qs::qread(paste0(dir_data,'/hic_df.qs'))


cat("Reshaping the data matrix for Singular Value Decomposition... \n")
mean_thres = 0
var_thres = 0
band_select = "all"

hic_df<-hic_df%>%mutate(diag=binB-binA)

n = length(unique(hic_df$cell))
setDT(hic_df)
summarized_hic = hic_df[, .(agg_m = mean(count), agg_v = var(count)),
                        by = .(chrom, binA, binB)]
summarized_hic[is.na(agg_v), "agg_v"] = 0#log(0.0001)

if (band_select == "all"){
  summarized_hic = summarized_hic %>%
    filter(binA - binB != 0, agg_m > quantile(agg_m, mean_thres),
           agg_v >= quantile(agg_v, var_thres)) %>% select(chrom, binA, binB)
}

cell_names = unique(hic_df$cell)
input_mat = matrix(0, nrow = length(cell_names), ncol = nrow(summarized_hic))

for (i in 1:length(cell_names)) {
  output_cell = summarized_hic
  output_cell$count = 0
  temp = hic_df[cell == cell_names[i], ]
  temp = temp %>% filter(diag > 0, chrom %in% summarized_hic$chrom,
                         binA %in% summarized_hic$binA,
                         binB %in% summarized_hic$binB)
  setDT(output_cell)
  setDT(temp)
  output_cell = output_cell[temp, `:=`(count, i.count), on = .(chrom,
                                                                     binA, binB)]
  input_mat[i, ] = output_cell$count
}

rm(temp)
rm(output_cell)


qs::qsave(summarized_hic,'summarized_hic.qs')




cat("Conducting initial SVD... \n")

svd_res=RSpectra::svds(input_mat,k = exploration_rank)

nn_sv=svd_res$d[svd_res$d>0]
(nn_sv ) %>% log %>% plot




pdf(file = paste0(dir_data,"/svd_plot.pdf"),width = 7,height = 7) 

# Step 2: Create the plot with R code
(nn_sv ) %>% log %>% plot(.,main="singular value plot of unfolded Hi-C",ylab="log(lambda)",xlab="Rank")
abline(v=10,col=2)
abline(v=20,col=3)
abline(v=30,col=4)
abline(v=40,col=5)
legend('topright', legend=c("Rank=10", "Rank=20", "Rank=30", "Rank=40"),
       col=c(2,3,4,5),lty = 1,cex=1)
dev.off()


cat("Please specify the rank value based on the singular value 'svd_plot.pdf' and hit enter (skipping will give rank=30) : ")
Rank <- as.numeric(readLines(con="stdin", 1))
cat(Rank, "\n")



#check rank truncation first
if(is.na(Rank)){Rank=30}
saveRDS(Rank,paste0(dir_data,'/Rank.rds'))
##############################################






 cell_type=unique(hic_df$cell)






svd_fitted=svd_res$u[,1:Rank]%*%diag(svd_res$d[1:Rank])%*%t(svd_res$v[,1:Rank])
rm(svd_res)
if(only_zero_entries==TRUE){input_mat[input_mat==0]=svd_fitted[input_mat==0];rm(svd_fitted)}
if(only_zero_entries!=TRUE){input_mat=svd_fitted;rm(svd_fitted)}

rownames(input_mat)=cell_type
rm(hic_df)

input_mat[input_mat<0]=0

if(debias==TRUE){
    
cat("Debiasing started... (If the memory does not allow debiasing, set debias=FALSE in the config_file_preprocess.R file and run Muscle again)\n")    
    
tmp_wide=data.frame(summarized_hic,t(input_mat))
rm(summarized_hic)
rm(input_mat)
tmp_long <- melt( tmp_wide, id.vars = c("chrom","binA","binB"))
rm(tmp_wide)
colnames(tmp_long)=c("chrom","binA","binB","cell","count")
tmp_df=tmp_long%>%mutate(diag=binB-binA)
rm(tmp_long)
hic_df_svd=tmp_df[,c(1:3,5,6,4)]
rm(tmp_df)
  #log transformation for debiasing
  hic_df_svd=hic_df_svd%>%mutate(count=log(count+0.0001))
  hic_df_svd = data.table(hic_df_svd)  
  band_info <-hic_df_svd[,.(band_depth = mean(count)),by = .(chrom,diag,cell)]
  alpha_j <- band_info[,.(depth = mean(band_depth)),by = .(chrom,diag)]
  hic_df_svd <- hic_df_svd %>% left_join(alpha_j, by = c("chrom", "diag")) 
  rm(alpha_j)    
    
hic_df_svd <- hic_df_svd %>% left_join(band_info,by = c("chrom", "diag", "cell")) %>% mutate(count = count-band_depth +depth) 
  rm(band_info)    
hic_df_svd <- hic_df_svd %>%select(-c(band_depth, depth))




n = length(unique(hic_df_svd$cell))

setDT(hic_df_svd)
summarized_hic = hic_df_svd[, .(agg_m = mean(count), agg_v = var(count)),
                            by = .(chrom, binA, binB)]
summarized_hic[is.na(agg_v), "agg_v"] =log(0.0001)


  summarized_hic = summarized_hic %>%
    filter(binA - binB != 0, agg_m > quantile(agg_m, mean_thres),
           agg_v >= quantile(agg_v, var_thres)) %>% select(chrom, binA, binB)

# quantile(agg_m, mean_thres) deletes locus pair attaining minimum (if non-transformed, delete zeros)

cell_names = unique(hic_df_svd$cell)
input_mat = matrix(log(0.0001), nrow = length(cell_names), ncol = nrow(summarized_hic))

for (i in 1:length(cell_names)) {
  output_cell = summarized_hic
  output_cell$count = log(0.0001)
  temp = hic_df_svd[cell == cell_names[i], ]
  temp = temp %>% filter(diag > 0, chrom %in% summarized_hic$chrom,
                         binA %in% summarized_hic$binA,
                         binB %in% summarized_hic$binB)
  setDT(output_cell)
  setDT(temp)
  output_cell = output_cell[temp, `:=`(count, i.count), on = .(chrom,
                                                                     binA, binB)]
  input_mat[i, ] = output_cell$count
}
rm(temp)
rm(hic_df_svd)

input_mat=exp(input_mat)

input_mat[input_mat==exp(log(0.0001))]=0

}


qs::qsave(summarized_hic,'summarized_hic.qs')

chrlist=paste0("chr",c(1:(chr_num-1),"X"))
for(chr in 1:chr_num){
  qs::qsave(input_mat[,which(summarized_hic$chrom==chrlist[chr])],paste0(dir_data,"/Matricized_HiCtensor_imputed_chr",chr,".qs"))
  print(chr)
  
}
rm(input_mat)
rm(summarized_hic)



print("Generating tensors...")

which.chrs=paste0(c(1:chr_num),collapse =  " ")
 
if(GNU==TRUE){
  if(!is.null(ssh))system(paste0('cd ',dir_functions,'; parallel -S ',ssh,' --jobs 4 --workdir . Rscript ::: tensor_generator.R ::: ', which.chrs,' ::: ',dir_data," ::: ",dir_functions," ::: ",chr_num," ::: ",sizefile))
  if(is.null(ssh))system(paste0('cd ',dir_functions,'; parallel --jobs 4 --workdir . Rscript ::: tensor_generator.R ::: ', which.chrs,' ::: ',dir_data," ::: ",dir_functions," ::: ",chr_num," ::: ",sizefile))
} 
if(GNU!=TRUE){
  for(chr in 1:chr_num){
    
    system(paste0('cd ',dir_functions,"; Rscript --vanilla tensor_generator.R ",chr,' ',' ',dir_data,' ',dir_functions,' ',chr_num,' ',sizefile))
  }
  
}
