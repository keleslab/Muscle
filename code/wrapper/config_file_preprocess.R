#number of chromosomes (23(Human) or 20(Mouse))
chr_num=20


#directory of the data
dir_data="{Muscle directory}/data/Liu2021"

#directory of the functions
dir_functions="{Muscle directory}/code/functions"


#Initial maximum rank of SVD of HiC
exploration_rank=150



#Ranks for each data 
#( 0 (default) : it means to choose rank based on the screeplot of SVD, e.g., Rank_HiC=0 means that the rank for Hi-C data is chosen based on SVD of matricized Hi-C.)
#( R : it means to set the rank as R, e.g., Rank_HiC=30 means that the rank for Hi-C data is 30. R is nonzero in this case.)

Rank_HiC=0
Rank_mCG=0
Rank_mCH=0



#Do you want debiasing? (Yes: TRUE, No: FALSE)
debias=TRUE

#Replace only zero entries with impuation (Yes: TRUE, No: FALSE)
only_zero_entries=FALSE


#If GNU parllel exists, type TRUE. Otherwise, FALSE. 
GNU=FALSE


#server names for gnu parallel. 
#If you are not using multiple servers for GNU parallel, leave it as ssh=NULL
ssh=NULL

#If multiple servers can be used, list them with comma without space. e.g., 'server01,server02'
#ssh='server01,server02,server03'


#chromosome size file name within the data directory.
#sizefile='mm9.chrom.sizes'
sizefile='mm10.chrom.sizes'
#sizefile='hg19.chrom.sizes'
