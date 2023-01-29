

#number of chromosomes (23(Human) or 20(Mouse))
chr_num=20


#directory of the data
dir_data="{Muscle directory}/data/example"

#directory of the functions
dir_functions="{Muscle directory}/code/functions"


#Initial maximum rank of SVD of HiC
exploration_rank=150

#Do you want debiasing? (Yes: TRUE, No: FALSE)
debias=TRUE

#Replace only zero entries with impuation (Yes: TRUE, No: FALSE)
only_zero_entries=FALSE


#If GNU parllel exists, type TRUE. Otherwise, FALSE. 
GNU=TRUE


#server names for gnu parallel. If exists, list them with comma without space. e.g., server01,server02
#ssh=NULL
ssh='hodor01,hodor02,hodor03'


#chromosome size file name within the data directory.
sizefile='mm10.chrom.sizes'
#sizefile='mm9.chrom.sizes'
#sizefile='hg19.chrom.sizes'

