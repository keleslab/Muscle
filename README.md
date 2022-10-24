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

For **config_file_preprocess.R**, only change **dir_data** and **dir_functions**. The [Muscle directory] is going to be the directory that you cloned the Muscle.e.g., /Users/kwangmoonpark/Muscle.

```
#number of chromosomes (23(Human) or 20(Mouse))
chr_num=20


#directory of the data
dir_data="[Muscle directory]/data/example"

#directory of the functions
dir_functions="[Muscle directory]/code/functions"


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

