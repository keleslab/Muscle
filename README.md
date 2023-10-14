# Muscle : semi-non negative joint decomposition of multiple single cell tensors
![Muscle diagram](/figures/Figure_intro.jpg)

## Muscle Usage

### 1. Preparation


First of all, before cloning the muscle github package, go to the right directory that you would like to implement Muscle. In the cmd terminal, do

```
cd {Muscle directory}
```

then go on to the next step. {Muscle directory} could be /Users/kp223/Muscle.


#### 1. Repository clone

For cloning the github repository, again on the cmd terminal, run the linux code 

```
git clone git@github.com:keleslab/Muscle.git
```

For R, we need the requirements as below : 


-   R: [R installation](https://www.r-project.org)  (>=4.2.1)
-   parallel (Not required, but highly recommended for faster implementation)

#### 2. Install/load required R packages

In R, run those codes that download the required packages for running Muscle.

```
install.packages('pacman')
pacman::p_load(MASS,Matrix,dplyr,rTensor,reshape2,Rcpp,RcppArmadillo,foreach,inline,parallel,doParallel,RSpectra,qs,gtools)
```

Details about implementing codes can be found in the [Wiki page](https://github.com/keleslab/Muscle/wiki) of this github. For the users with low memory or the users who wants to just check if Muscle works, we recommend to first try this [Li et al 2019 data tutorial](https://github.com/keleslab/Muscle/wiki/Li-et-al-2019).

In case you are already done with installing Muscle, you can direclty move onto section 2 of each of the Wiki pages.
