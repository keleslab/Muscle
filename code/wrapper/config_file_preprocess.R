

#number of chromosomes (23(Human) or 20(Mouse))
chr_num=20


#directory of the data
dir_data="/storage08/kwangmoon/Muscle/data/example"

#directory of the functions
dir_functions="/storage08/kwangmoon/Muscle/code/functions"


#Initial maximum rank of SVD
exploration_rank=150

#debias?
debias=FALSE

#Replace only zero entries with impuation
only_zero_entries=TRUE

#server names for gnu parallel. If exists, list them with comma without space. e.g., server01,server02
ssh=NULL

#chromosome size file name within the data directory.
sizefile='mm9.chrom.sizes'
#sizefile='hg19.chrom.sizes'

