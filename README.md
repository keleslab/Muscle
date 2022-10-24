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
This section gives instruction for running the example [Li \textit{et al.} 2019 data](https://www.nature.com/articles/s41592-019-0502-z) 
First change directory to **wrapper** as below.

```
cd /Users/kwangmoonpark/Muscle/code/wrapper/
```
Before running **Muscle.sh** file, change two configuration files in **wrapper** directory: **config_file_preprocess.R**, **config_file_model.R**, which contain aruments and hyperparameters for running Muscle.

For **config_file_preprocess.R**, 


    1. chr_num : 
