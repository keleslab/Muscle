# 
# chr_num=20
# dir='/storage08/kwangmoon'
# exploration_rank=100
# debias=FALSE
# only_zero_entries=TRUE
# ssh=NULL
# sizefile='mm9.chrom.sizes'
# 
# 

# 
# 
# args = commandArgs(trailingOnly=TRUE)
# 
# ######
# 
# chr_num=as.numeric(args[1])
# dir=(args[2])
# exploration_rank=as.numeric(args[3])
# debias=as.logical(args[4])
# only_zero_entries=as.logical(args[5])
# ssh=as.character(args[6])
# sizefile=as.character(args[7])


source("config_file_preprocess.R")


#######
# dir_data=paste0(dir,"/data/example")
# dir_functions=paste0(dir,"/code/functions")
setwd(dir_data)

saveRDS(chr_num,paste0(dir_data,'/chr_num.rds'))

pacman::p_load(RSpectra,qs,reshape2,dplyr,rTensor,Matrix,data.table)
options(scipen = 4)


##################################3

hic_df=qs::qread(paste0(dir_data,'/hic_df.qs'))

mean_thres = 0
var_thres = 0
band_select = "all"



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
print(paste("The number of features is", nrow(summarized_hic)))
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




# qs::qsave(input_mat,'longform_hic.qs')
qs::qsave(summarized_hic,'summarized_hic.qs')
# 
# 
# input_mat=qs::qread('longform_hic.qs')
# summarized_hic=qs::qread('summarized_hic.qs')



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



#Rank=readline("Please specify the rank value based on the singular value plot pdf and hit enter (skipping will give rank=30) : ");print(Rank)


#qs::qsave(svd_res,"svd_res_raw_Li2019_rankfull.qs")
#check rank truncation first
if(is.na(Rank)){Rank=30}
saveRDS(Rank,paste0(dir_data,'/Rank.rds'))
##############################################





# hic_df=qs::qread("Li2019.qs")
 cell_type=unique(hic_df$cell)#qs::qread('cell_type.qs')


# 
# summarized_hic=qs::qread('summarized_hic_Li2019.qs')
# svd_res=qs::qread("svd_res_raw_Li2019_rankfull.qs")


svd_fitted=svd_res$u[,1:Rank]%*%diag(svd_res$d[1:Rank])%*%t(svd_res$v[,1:Rank])

if(only_zero_entries==TRUE){input_mat[input_mat==0]=svd_fitted[input_mat==0];rm(svd_fitted)}
if(only_zero_entries!=TRUE){input_mat=svd_fitted;rm(svd_fitted)}
#cell_type=unique(hic_df$cell)#qs::qread('cell_type.qs')
rownames(input_mat)=unique(hic_df$cell)


input_mat[input_mat<0]=0

if(debias==TRUE){
tmp_wide=data.frame(summarized_hic,t(input_mat))
tmp_long <- melt( tmp_wide, id.vars = c("chrom","binA","binB"))
colnames(tmp_long)=c("chrom","binA","binB","cell","count")
tmp_df=tmp_long%>%mutate(diag=binB-binA)

hic_df_svd=tmp_df[,c(1:3,5,6,4)]

  #log transformation for debiasing
  hic_df_svd=hic_df_svd%>%mutate(count=log(count+0.0001))
  band_info <- hic_df_svd%>% group_by(chrom, diag, cell) %>% summarise(band_depth = mean(count))
  alpha_j <- band_info %>% group_by(chrom, diag) %>% summarise(depth = mean(band_depth))
  hic_df_svd <- hic_df_svd %>% left_join(alpha_j, by = c("chrom", "diag")) %>% left_join(band_info,
                                                                                         by = c("chrom", "diag", "cell")) %>% mutate(count = count-band_depth +depth) %>% 
    select(-c(band_depth, depth, count))
  
  





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
print(paste("The number of features is", nrow(summarized_hic)))
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



input_mat=exp(input_mat)

input_mat[input_mat==exp(log(0.0001))]=0

}


qs::qsave(summarized_hic,'summarized_hic.qs')

chrlist=paste0("chr",c(1:(chr_num-1),"X"))
for(chr in 1:chr_num){
  qs::qsave(input_mat[,which(summarized_hic$chrom==chrlist[chr])],paste0(dir_data,"/Matricized_HiCtensor_imputed_chr",chr,".qs"))
  print(chr)
  
}



print("Generating tensors...")

which.chrs=paste0(c(1:chr_num),collapse =  " ")
 
 
if(!is.null(ssh))system(paste0('cd ',dir_functions,'; parallel -S ',ssh,' --jobs 4 --workdir . Rscript ::: tensor_generator.R ::: ', which.chrs,' ::: ',dir_data," ::: ",dir_functions," ::: ",chr_num," ::: ",sizefile))
if(is.null(ssh))system(paste0('cd ',dir_functions,'; parallel --jobs 4 --workdir . Rscript ::: tensor_generator.R ::: ', which.chrs,' ::: ',dir_data," ::: ",dir_functions," ::: ",chr_num," ::: ",sizefile))

