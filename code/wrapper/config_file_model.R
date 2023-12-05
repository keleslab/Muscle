


#Directory of data: The directory should contain hic_df.qs and chrom.sizes file.
#dir_data="/storage08/kwangmoon/Muscle/data/Tan2021"
dir_data="/storage08/kwangmoon/Muscle/data/Lee2019/500kb"

#Directory of functions
dir_functions="/storage08/kwangmoon/Muscle/code/functions"

#Directory where Muscle output goes into

#dir_out="/storage08/kwangmoon/Muscle/results/Tan2021"
dir_out="/storage08/kwangmoon/Muscle/results/Lee2019/500kb/schic"



#Initial maximum rank of Methylation SVD
exploration_rank=150


#Modality of the tensors. #(only 3 choices available : "All", or "HiC", or "HiC+CG")

#modality="All" #When all Hi-C, mCG, mCH are analyzed
modality="HiC" #When only Hi-C is analyzed
#modality="HiC+CG" #when only Hi-C and mCG are analyzed



#Ranks for each data 
#( 0 (default) : it means to choose rank based on the screeplot of SVD, e.g., Rank_HiC=0 means that the rank for Hi-C data is chosen based on SVD of Hi-C.)
#( R : it means to set the rank as R, e.g., Rank_HiC=30 means that the rank for Hi-C data is 30. R is nonzero in this case.)

Rank_HiC=50
Rank_mCG=30
Rank_mCH=0




# If bulk TAD information exists, set it as TRUE. Otherwise, FALSE
#Bulk_exist=TRUE
Bulk_exist=FALSE

#If GNU parllel exists, type TRUE. Otherwise, FALSE. 
GNU=TRUE
#GNU=FALSE

#tolerance level for each rank 1 update 
tol=0.00001


#maximum number of iterations for each rank 1 update
maxiter=10

#server names for gnu parallel. 
#If you are not using multiple servers for GNU parallel, leave it as ssh=NULL
#If multiple servers can be used, list them with comma without space. 
#e.g., ssh='server01,server02,server03'
ssh='hodor01,hodor02'
