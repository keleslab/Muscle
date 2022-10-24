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

First change directory to wrapper as below.

```
cd /Users/kwangmoonpark/Muscle/code/wrapper/
```
Before running **Muscle.sh** file, 

