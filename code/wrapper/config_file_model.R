


#Directory of data: The directory should contain hic_df.qs and chrom.sizes file.
dir_data="{Muscle directory}/data/Liu2021"

#Directory of functions
dir_functions="{Muscle directory}/code/functions"

#Directory where Muscle output goes into

dir_out="{Muscle directory}/results/Liu2021"



#Initial maximum rank of Methylation SVD
exploration_rank=150


#Modality of the tensors. #(only 3 choices available : "All", or "HiC", or "HiC+CG")

modality="All" #When all Hi-C, mCG, mCH are analyzed
#modality="HiC" #When only Hi-C is analyzed
#modality="HiC+CG" #when only Hi-C and mCG are analyzed


# If bulk TAD information exists, set it as TRUE. Otherwise, FALSE
#Bulk_exist=TRUE
Bulk_exist=FALSE

#If GNU parllel exists, type TRUE. Otherwise, FALSE. 
#GNU=TRUE
GNU=FALSE

#tolerance level for each rank 1 update 
tol=0.00001


#maximum number of iterations for each rank 1 update
maxiter=10

#server names for gnu parallel. 
#If you are not using multiple servers for GNU parallel, leave it as ssh=NULL
ssh=NULL

#If multiple servers can be used, list them with comma without space. e.g., 'server01,server02'
#ssh='server01,server02,server03'
