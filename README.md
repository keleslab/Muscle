# Muscle : semi-non negative joint decomposition of multiple single cell tensors
![Muscle diagram](/figures/Figure_intro.jpg)

## Muscle Usage

### 1. Preparation

git clone git@github.com:keleslab/Muscle.git

-   R: [R installation](https://www.r-project.org)  (>=4.2.1)
-   parallel (Not required, but highly recommended)

#### 1. Required R packages

```
install.packages('pacman')
pacman::p_load(MASS,Matrix,dplyr,rTensor,reshape2,Rcpp,foreach,inline,parallel,doParallel,RSpectra,qs)
```


#### 2. Quick start
This section gives instruction for running the example [Li et al. 2019 data](https://www.nature.com/articles/s41592-019-0502-z) 



First change directory to **wrapper** as below.

```
cd /Users/kwangmoonpark/Muscle/code/wrapper/
```
Before running **Muscle.sh** file, change two configuration files in **wrapper** directory: **config_file_preprocess.R**, **config_file_model.R**, which contain aruments and hyperparameters for running Muscle.

For **config_file_preprocess.R**, only change **dir_data** and **dir_functions**. The {Muscle directory} is going to be the directory that you cloned the Muscle. e.g., {Muscle directory} can be /Users/kwangmoonpark/Muscle.

In case **parallel** does not exist, one can set **GNU=FALSE**. Moreover, if multiple servers are available for GNU parallel, one can set ssh argument such as **ssh='server01,server02,server03'**. Note that **ssh=NULL** should be used when **parallel** does not exist.


```
#number of chromosomes (23(Human) or 20(Mouse))
chr_num=20


#directory of the data
dir_data="{Muscle directory}/data/example"

#directory of the functions
dir_functions="{Muscle directory}/code/functions"


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

```




For **config_file_model.R**, only change **dir_data**,**dir_functions**, **dir_out**. The {Muscle directory} is going to be the directory that you cloned the Muscle. e.g., {Muscle directory} can be /Users/kwangmoonpark/Muscle.

Similar to **config_file_preprocess.R** file, in case **parallel** does not exist, one can set **GNU=FALSE**. Moreover, if multiple servers are available for GNU parallel, one can set ssh argument such as **ssh='server01,server02,server03'**. Note that **ssh=NULL** should be used when **parallel** does not exist.


```



#Directory of data: The directory should contain hic_df.qs and chrom.sizes file.
dir_data="{Muscle directory}/data/example"

#Directory of functions
dir_functions="{Muscle directory}/code/functions"

#Directory where Muscle output goes into
dir_out="{Muscle directory}/results/Li2019"



#Initial maximum rank of Methylation SVD
exploration_rank=150


#Modality of the tensors. #Only threecasesare allowed (All,HiC,HiC+CG)

#modality="All" #When all Hi-C, mCG, mCH are analyzed
#modality="HiC" #When only Hi-C is analyzed
modality="HiC+CG" #when only Hi-C and mCG are analyzed


# If bulk TAD information exists, set it as TRUE. Otherwise, FALSE
#Bulk_exist=TRUE
Bulk_exist=FALSE

#If GNU parllel exists, type TRUE. Otherwise, FALSE. 
GNU=TRUE

#tolerance level for each rank 1 update 
tol=0.00001


#maximum number of iterations for each rank 1 update
maxiter=4

#servers should be listed without any space in between. If you are not using multiple servers for GNU parallel, leave it as ssh=NULL
#For example, ssh="hodor01,hodor02,hodor03"
ssh=NULL


```


