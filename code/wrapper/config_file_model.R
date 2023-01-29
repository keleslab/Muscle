


#Directory of data: The directory should contain hic_df.qs and chrom.sizes file.
dir_data="{Muscle directory}/data/example"

#Directory of functions
dir_functions="{Muscle directory}/code/functions"

#Directory where Muscle output goes into

dir_out="{Muscle directory}/results/example"



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
GNU=TRUE

#tolerance level for each rank 1 update 
tol=0.00001


#maximum number of iterations for each rank 1 update
maxiter=10

#servers should be listed without any space in between. If you are not using multiple servers for GNU parallel, leave it as ssh=NULL
#For example, ssh="hodor01,hodor02,hodor03"
ssh="hodor01,hodor02,hodor03"




