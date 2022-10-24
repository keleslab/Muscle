



#Directory of data: The directory should contain hic_df.rds and chrom.sizes file.
dir_data="/storage08/kwangmoon/Muscle/data/example"

#Directory of functions
dir_functions="/storage08/kwangmoon/Muscle/code/functions"

#Directory where Muscle output goes into
dir_out="/storage08/kwangmoon/Muscle/results/Li2019"



#Initial maximum rank of SVD
exploration_rank=150


#Rank of HiC tensors
Rank=as.numeric(readRDS(paste0(dir_data,'/Rank.rds')))
 

#Modality of the tensors. 
#Only threecasesare allowed (All,HiC,HiC+CG)

#modality="All" #When all Hi-C, mCG, mCH are analyzed
# modality="HiC" #When only Hi-C is analyzed
modality="HiC+CG" #when only Hi-C and mCG are analyzed


# number of chromosomes(human:23, mouse:20)
chr_num=as.numeric(readRDS(paste0(dir_data,'/chr_num.rds')))



# If bulk TAD information exists, set it as TRUE. Otherwise, FALSE
#Bulk_exist=TRUE
Bulk_exist=FALSE



#tolerance level for each rank 1 update 
tol=0.00001


#maximum number of iterations for each rank 1 update
maxiter=4

#servers should be listed without any space in between. If you are not using multiple servers for GNU parallel, leave it as ssh=NULL
#For example, ssh="hodor01,hodor02,hodor03"
ssh=NULL

