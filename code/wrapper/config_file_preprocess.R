

#number of chromosomes (23(Human) or 20(Mouse))
chr_num=20


#directory of the data
dir_data="/Users/kwangmoonpark/Muscle/data/example"

#directory of the functions
dir_functions="/Users/kwangmoonpark/Muscle/code/functions"


#Initial maximum rank of SVD of HiC
exploration_rank=150

#debias?
debias=FALSE

#Replace only zero entries with impuation
only_zero_entries=TRUE

#server names for gnu parallel. If exists, list them with comma without space. e.g., server01,server02
ssh=NULL

#If GNU parllel exists, type TRUE. Otherwise, FALSE. 
GNU=TRUE


#chromosome size file name within the data directory.
sizefile='mm9.chrom.sizes'
#sizefile='hg19.chrom.sizes'

