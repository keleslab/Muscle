###################
##Matrix operators#
###################


BTD_combi2=function(A1,A2,A3,R,Lr){
  Theta=0
  for (r in 1:R) {
    
    if(r==1)Theta=(A1[,1:Lr[r]]%*%t(A2[,1:Lr[r]]))%o%A3[,r] 
    else{
      
      Theta=Theta+(A1[,(sum(Lr[1:(r-1)])+1):sum(Lr[1:r])]%*%t(A2[,(sum(Lr[1:(r-1)])+1):sum(Lr[1:r])]))%o%A3[,r]   
      
    }
    
  }
  return(Theta)
}


khatri_rao_matrix=function(a,b,R,la,lb){
  
  #  blksize_a=dim(a)[2]/R  
  #  blksize_b=dim(b)[2]/R  
  
  res=c()
  for(r in 1:R){
    
    #res=cbind(res,kronecker(a[,(blksize_a*r-(blksize_a-1)):(blksize_a*r)],b[,(blksize_b*r-(blksize_b-1)):(blksize_b*r)]))
    if(r==1) res=cbind(res,kronecker(a[,1:la[r]],b[,1:lb[r]]))
    else{
      
      res=cbind(res,kronecker(a[,(sum(la[1:(r-1)])+1):sum(la[1:r])],b[,(sum(lb[1:(r-1)])+1):sum(lb[1:r])]))
      
    }
    
  }
  
  
  return(res)
}



khatri_rao_matrix2=function(a,b,R,Lr){
  #  blksize_a=dim(a)[2]/R  
  #  blksize_b=dim(b)[2]/R  
  
  res=c()
  for(r in 1:R){
    
    # res=cbind(res,khatri_rao(a[,(blksize_a*r-(blksize_a-1)):(blksize_a*r)],b[,(blksize_b*r-(blksize_b-1)):(blksize_b*r)])%*%rep(1,Lr[r]))  
    if(r==1)res=cbind(res,khatri_rao(a[,1:Lr[r]],b[,1:Lr[r]])%*%rep(1,Lr[r]))  
    else{
      
      res=cbind(res,khatri_rao(a[,(sum(Lr[1:(r-1)])+1):sum(Lr[1:r])],b[,(sum(Lr[1:(r-1)])+1):sum(Lr[1:r])])%*%rep(1,Lr[r]))  
      
    }  
  }
  
  
  return(res)
}




Fnorm=function(x){sqrt(sum((x)^2))}

##############
#data process#
##############

naomit_tensor=function(x){
  
  
  return(x[which(rowMeans(x!=0)!=0), which(colMeans(x!=0)!=0)])
  
}
naomit_matrix=function(x){
  
  
  return(x[which(rowMeans(x!=0)!=0), ])
  
}

scalemat=function(x){(x-mean(x))/sd(x)}


####
#tensor generators
####



sparse2dense <- function(sparse.as.df, N = NULL, balancing = TRUE, missing.cells.na = FALSE){
  colnames(sparse.as.df) <- c("i","j","val")
  with(
    sparse.as.df,
    Matrix::sparseMatrix(i=as.numeric(i),
                         j=as.numeric(j),
                         x=as.numeric(val))
  ) %>%
    as.matrix() -> dense.mtx
  if((!is.null(N)) | balancing){
    dense.mtx <- balance(dense.mtx, N = N)
  }
  shape <- dim(dense.mtx)
  if(
    (shape[1] == shape[2]) &
    (all(sparse.as.df$i >= sparse.as.df$j) | all(sparse.as.df$i <= sparse.as.df$j))
  ){
    # if matrix is square symmetric then symmetrize it
    dense.mtx.diag <- diag(dense.mtx)
    dense.mtx <- dense.mtx + t(dense.mtx)
    diag(dense.mtx) <- dense.mtx.diag
  }
  # lastly change 0 to NA if missing.cells.na
  if(missing.cells.na & any(dense.mtx == 0)){
    # convert missing cells to NA instead of 0
    if.missing <- matrix(0, nrow = nrow(dense.mtx), ncol = ncol(dense.mtx))
    if.missing <- if.missing == 0
    # set non missing cells to FALSE
    if.missing[cbind(sparse.as.df$i, sparse.as.df$j)] <- FALSE
    dense.mtx[if.missing] <- NA
  }
  return(dense.mtx)
}

balance <- function(mtx, N = NULL){
  if(!is.null(N)){
    stopifnot(N > 0)
  }
  mtx.n.rows <- nrow(mtx)
  mtx.n.cols <- ncol(mtx)
  # balance if nrows != ncols
  if(mtx.n.rows < mtx.n.cols){
    # add 0-s rows at the end
    rep(0, mtx.n.cols) %>%
      replicate(mtx.n.cols - mtx.n.rows,.) %>%
      t() %>% rbind(mtx, .) -> mtx
  } else if(mtx.n.rows > mtx.n.cols){
    # add 0-s columns at the end
    rep(0, mtx.n.rows) %>%
      replicate(mtx.n.rows - mtx.n.cols,.) %>%
      cbind(mtx, .) -> mtx
  }
  # add/remove columns and rows, so dim(mtx) == c(N,N)
  if(!is.null(N)){
    if(N > nrow(mtx)){
      rep(0, nrow(mtx)) %>%
        replicate(N - nrow(mtx), .) %>%
        cbind(mtx, .) -> mtx
      rep(0, ncol(mtx)) %>%
        replicate(N - nrow(mtx), .) %>%
        t() %>% rbind(mtx, .) -> mtx
    } else if(N < nrow(mtx)){
      mtx <- mtx[1:N,1:N]
    }
  }
  return(mtx)
}


######
# Initial hosvd Tensor decomp
######



hosvd_new <- function(tnsr,ranks=NULL){
  
  num_modes <- tnsr@num_modes
  #no truncation if ranks not provided
  #progress bar
  pb <- txtProgressBar(min=0,max=num_modes,style=3)
  #loops through and performs SVD on mode-m matricization of tnsr
  U_list <- vector("list",num_modes)
  for(m in 1:num_modes){
    temp_mat <- rTensor::rs_unfold(tnsr,m=m)@data
    U_list[[m]] <- RSpectra::svds(temp_mat,k=ranks[m])$u
    setTxtProgressBar(pb,m)
  }
  close(pb)
  #computes the core tensor
  Z <- rTensor::ttl(tnsr,lapply(U_list,t),ms=1:num_modes)
  
  #put together the return list, and returns
  list(Z=Z,U=U_list)
}

######
#JBTD#
######
rankone_Muscle=function(data,k,LR,tol=0.0001,maxiter=50,dir_out,dir_functions,chr_num,modality,ssh="NULL"){
  R=1
  Lr_list=as.list(LR)
  
  dim_list=lapply(data, dim)
  d_list=lapply(dim_list, prod)
  Nten=d_list[1:chr_num] %>% unlist() %>% sum
  if(modality=="All"){d_list=as.list(c(rep(Nten,chr_num),d_list[[chr_num+1]],d_list[[chr_num+2]]))}
  if(modality=="HiC+CG"){d_list=as.list(c(rep(Nten,chr_num),d_list[[chr_num+1]]))}
  if(modality=="HiC"){d_list=as.list(c(rep(Nten,chr_num)))}
  
  
  
  ######
  btd_res=list()
  for(chr in 1:chr_num){
    btd_res[[chr]]=list(A=c(),B=c(),C=c())
  }
  
  if(modality=="All"|modality=="HiC+CG"){
    methyl_res_CG=list(Methyl_loci_loading=c(),Methyl_cell_loading=c())
    if(modality=="All"){methyl_res_CH=list(Methyl_loci_loading=c(),Methyl_cell_loading=c())}
    
    
  }
  
  
  cat("Initializing...")
  cat("\n")
  
  #  which.cells=paste0(c(1:chr_num),collapse =  " ")
  if(modality=="All"){which.tensors=paste0(c(1:chr_num,"CG","CH"),collapse =  " ")}
  if(modality=="HiC"){which.tensors=paste0(c(1:chr_num),collapse =  " ")}
  if(modality=="HiC+CG"){which.tensors=paste0(c(1:chr_num,"CG"),collapse =  " ")}
  
  
  outname=paste0(dir_out,'/data')
  
  
  ########################
  
  
  if(!is.null(ssh))system(paste0('cd ',dir_functions,'; parallel -S ',ssh,' --jobs 4 --workdir . Rscript ::: Initializer.R ::: ', which.tensors, ' ::: ',dir_out))
  if(is.null(ssh))system(paste0('cd ',dir_functions,'; parallel --jobs 4 --workdir . Rscript ::: Initializer.R ::: ', which.tensors, ' ::: ',dir_out))
  
  
  
  for(chr in 1:chr_num){
    
    
    btd_res[[chr]]$A=tryCatch(readRDS(paste0(outname,chr,'_modeA_res.rds')),error=function(e){system(paste0('cd ',dir_functions,"; Rscript --vanilla Initializer.R ",
                                                                                                            chr,' ', dir_out));return(readRDS(paste0(outname,chr,'_modeA_res.rds')))})
    btd_res[[chr]]$B=readRDS(paste0(outname,chr,'_modeB_res.rds'))
    
  }
  
  
  if(modality=="All"|modality=="HiC+CG"){
           chr="CG"
            Methyl_loci_loading_CG=tryCatch(readRDS(paste0(outname,chr,'_Methyl_loci_res.rds')),error=function(e){system(paste0('cd ',dir_functions,"; Rscript --vanilla Initializer.R ",
                                                                                                                        chr,' ', dir_out));return(readRDS(paste0(outname,chr,'_Methyl_loci_res.rds')))})
      if(modality=="All"){
            chr="CH"
          Methyl_loci_loading_CH=tryCatch(readRDS(paste0(outname,chr,'_Methyl_loci_res.rds')),error=function(e){system(paste0('cd ',dir_functions,"; Rscript --vanilla Initializer.R ",
                                                                                                                      chr,' ', dir_out));return(readRDS(paste0(outname,chr,'_Methyl_loci_res.rds')))})
                      }    
  
    
    
                    }
  
  
  junk <- dir(path=dir_out,  pattern="_mode")
  file.remove(paste0(dir_out,"/",junk))
  if(modality=="All"|modality=="HiC+CG"){
    junk <- dir(path=dir_out,  pattern="_Methyl_loci_res")
    file.remove(paste0(dir_out,"/",junk))  
  }
  
  
  #######3
  
  #k=3
  
  
  if(modality=="All"){btd_res[[1]]$C=t(as.matrix(Reduce('+',mapply(function(y,d){1/d*readRDS(paste0(outname,y,'_XTY.rds'))},as.list(c(1:chr_num,"CG","CH")),d_list))))}
  if(modality=="HiC+CG"){btd_res[[1]]$C=t(as.matrix(Reduce('+',mapply(function(y,d){1/d*readRDS(paste0(outname,y,'_XTY.rds'))},as.list(c(1:chr_num,"CG")),d_list))))}
  if(modality=="HiC"){btd_res[[1]]$C=t(as.matrix(Reduce('+',mapply(function(y,d){1/d*readRDS(paste0(outname,y,'_XTY.rds'))},as.list(c(1:chr_num)),d_list))))}
  
  
  
  embeddings=btd_res[[1]]$C
  
  
  mult=1
  if(all(embeddings[,1]<0)){mult=-1;print(paste0(1,"Module sign flipped!"))}
  
  btd_res[[1]]$C[,1]=mult*btd_res[[1]]$C[,1]
  btd_res[[1]]$C[,1]=ifelse(btd_res[[1]]$C[,1]<0,0,btd_res[[1]]$C[,1])
  
  for(chr in 1:chr_num){
    ind=(Lr_list[[chr]][1]*1-(Lr_list[[chr]][1]-1)):(Lr_list[[chr]][1]*1) 
    
    btd_res[[chr]]$A[,ind]=mult*btd_res[[chr]]$A[,ind]
    
    
  }      
  
  if(modality=="All"|modality=="HiC+CG"){
    Methyl_loci_loading_CG=mult*Methyl_loci_loading_CG
    if(modality=="All"){Methyl_loci_loading_CH=mult*Methyl_loci_loading_CH}
    
  }
  
  
  for(chr in 1:chr_num){
    
    btd_res[[chr]]$C=btd_res[[1]]$C
    
  }
  
  if(modality=="All"|modality=="HiC+CG"){
  Methyl_cell_loading_CG=btd_res[[1]]$C
  if(modality=="All"){Methyl_cell_loading_CH=btd_res[[1]]$C}
  }
  
  for( chr in 1:chr_num){btd_res[[chr]]$C=apply(btd_res[[chr]]$C,2,function(x)ifelse(is.finite(x/Fnorm(x)),x/Fnorm(x),0))}
  
  
  if(modality=="All"|modality=="HiC+CG"){
  Methyl_cell_loading_CG=apply(Methyl_cell_loading_CG,2,function(x)ifelse(is.finite(x/Fnorm(x)),x/Fnorm(x),0))
  if(modality=="All"){Methyl_cell_loading_CH=apply(Methyl_cell_loading_CH,2,function(x)ifelse(is.finite(x/Fnorm(x)),x/Fnorm(x),0))}
  }
  
  if(modality=="All"|modality=="HiC+CG"){
  methyl_res_CG$Methyl_loci_loading=Methyl_loci_loading_CG
  methyl_res_CG$Methyl_cell_loading=Methyl_cell_loading_CG
  
  if(modality=="All"){  methyl_res_CH$Methyl_loci_loading=Methyl_loci_loading_CH
  methyl_res_CH$Methyl_cell_loading=Methyl_cell_loading_CH}

  }
  Theta_list_new=mapply(function(x,y)BTD_combi2(x$A,x$B,x$C,R,y),btd_res,Lr_list,SIMPLIFY = F)
  
  if(modality=="All"|modality=="HiC+CG"){
  Theta_list_new[['CG']]=methyl_res_CG$Methyl_loci_loading%*%t(methyl_res_CG$Methyl_cell_loading)
  if(modality=="All"){Theta_list_new[['CH']]=methyl_res_CH$Methyl_loci_loading%*%t(methyl_res_CH$Methyl_cell_loading)}
  }
  ls_list=mapply(function(x,y){(sum((x-y)^2))},data,Theta_list_new,SIMPLIFY = F)

  
  
  
  f_new=Reduce('+',mapply(function(s,d){1/d*s},ls_list,d_list,SIMPLIFY = F))
  
  
  
  
  
  
  
  
  
  iter=0
  repeat {
    print(iter)
    
    
    

    
    outname=paste0(dir_out,'/data')
    
    
    
    
    cat('Saving iteration result...')
    saveRDS(btd_res,paste0(dir_out,'/HiC_result_iter_',iter,'_rank',k,'.rds'))
    if(modality=="All"|modality=="HiC+CG"){
    saveRDS(methyl_res_CG,paste0(dir_out,'/methyl_CG_result_iter_',iter,'_rank',k,'.rds'))
      if(modality=="All"){saveRDS(methyl_res_CH,paste0(dir_out,'/methyl_CH_result_iter_',iter,'_rank',k,'.rds'))      }
    }
    iter=iter+1
    
    ###k=1...3
    
    
    #k=1            
    
    
    cat('Estimating Theta...\n') 
    
    
    #k=1            
    
    btd_res_old=btd_res
    if(modality=="All"|modality=="HiC+CG"){
    methyl_res_CG_old=methyl_res_CG
    if(modality=="All"){methyl_res_CH_old=methyl_res_CH}
    }
    
    for(chr in 1:chr_num){
      
      saveRDS(btd_res[[chr]],paste0(outname,chr,'_btd_res.rds'))
      
    }
    
    if(modality=="All"|modality=="HiC+CG"){
    saveRDS(methyl_res_CG,paste0(outname,'CG','_Methyl_res.rds'))
      if(modality=="All"){saveRDS(methyl_res_CH,paste0(outname,'CH','_Methyl_res.rds'))}
    }
    
    #update mode A and B(symmetric)    
    if(!is.null(ssh))system(paste0('cd ',dir_functions,'; parallel -S ',ssh,' --jobs 4 --workdir . Rscript ::: ModeAB_learner.R ::: ', which.tensors, ' ::: ',dir_out))
    if(is.null(ssh))system(paste0('cd ',dir_functions,'; parallel --jobs 4 --workdir . Rscript ::: ModeAB_learner.R ::: ', which.tensors, ' ::: ',dir_out))
    
    
    for(chr in 1:chr_num){
      
      
      btd_res[[chr]]$A=tryCatch(readRDS(paste0(outname,chr,'_modeA_res.rds')),error=function(e){system(paste0('cd ',dir_functions,"; Rscript --vanilla ModeAB_learner.R ",
                                                                                                              chr,' ',dir_out));return(readRDS(paste0(outname,chr,'_modeA_res.rds')))})
      btd_res[[chr]]$B=readRDS(paste0(outname,chr,'_modeB_res.rds'))
      
    }
    
    if(modality=="All"|modality=="HiC+CG"){
    chr="CG"
    Methyl_loci_loading_CG=tryCatch(readRDS(paste0(outname,chr,'_Methyl_loci_res.rds')),error=function(e){system(paste0('cd ',dir_functions,"; Rscript --vanilla ModeAB_learner.R ",
                                                                                                                        chr,' ',dir_out));return(readRDS(paste0(outname,chr,'_Methyl_loci_res.rds')))})
    
    if(modality=="All"){    chr="CH"
    Methyl_loci_loading_CH=tryCatch(readRDS(paste0(outname,chr,'_Methyl_loci_res.rds')),error=function(e){system(paste0('cd ',dir_functions,"; Rscript --vanilla ModeAB_learner.R ",
                                                                                                                        chr,' ',dir_out));return(readRDS(paste0(outname,chr,'_Methyl_loci_res.rds')))})}

    
    
    
    }
    #files removing
    junk <- dir(path=dir_out,  pattern="_mode")
    file.remove(paste0(dir_out,"/",junk))
    if(modality=="All"|modality=="HiC+CG"){
      junk <- dir(path=dir_out,  pattern="_Methyl_loci_res")
      file.remove(paste0(dir_out,"/",junk))  
    }
    #3rd mode opt for initialization
    #######3
    
    #k=3
    
    
    if(modality=="All"){btd_res[[1]]$C=t(as.matrix(Reduce('+',mapply(function(y,d){1/d*readRDS(paste0(outname,y,'_XTY.rds'))},as.list(c(1:chr_num,"CG","CH")),d_list))))}
    if(modality=="HiC+CG"){btd_res[[1]]$C=t(as.matrix(Reduce('+',mapply(function(y,d){1/d*readRDS(paste0(outname,y,'_XTY.rds'))},as.list(c(1:chr_num,"CG")),d_list))))}
    if(modality=="HiC"){btd_res[[1]]$C=t(as.matrix(Reduce('+',mapply(function(y,d){1/d*readRDS(paste0(outname,y,'_XTY.rds'))},as.list(c(1:chr_num)),d_list))))}
    
    
    embeddings=btd_res[[1]]$C
    
    
    mult=1
    
    if(all(embeddings[,1]<0)){mult=-1;print(paste0(1,"Module sign flipped!"))}
    
    btd_res[[1]]$C[,1]=mult*btd_res[[1]]$C[,1]
    btd_res[[1]]$C[,1]=ifelse(btd_res[[1]]$C[,1]<0,0,btd_res[[1]]$C[,1])
    
    for(chr in 1:chr_num){
      ind=(Lr_list[[chr]][1]*1-(Lr_list[[chr]][1]-1)):(Lr_list[[chr]][1]*1) 
      
      btd_res[[chr]]$A[,ind]=mult*btd_res[[chr]]$A[,ind]
      
      
    }      
    if(modality=="All"|modality=="HiC+CG"){
    Methyl_loci_loading_CG=mult*Methyl_loci_loading_CG
    if(modality=="All"){Methyl_loci_loading_CH=mult*Methyl_loci_loading_CH}
    }
    
    for(chr in 1:chr_num){
      
      btd_res[[chr]]$C=btd_res[[1]]$C
      
    }
    if(modality=="All"|modality=="HiC+CG"){
      Methyl_cell_loading_CG=btd_res[[1]]$C
      if(modality=="All"){Methyl_cell_loading_CH=btd_res[[1]]$C  }
      
    }
    
    
    
    for( chr in 1:chr_num){btd_res[[chr]]$C=apply(btd_res[[chr]]$C,2,function(x)ifelse(is.finite(x/Fnorm(x)),x/Fnorm(x),0))}
    
    if(modality=="All"|modality=="HiC+CG"){
    Methyl_cell_loading_CG=apply(Methyl_cell_loading_CG,2,function(x)ifelse(is.finite(x/Fnorm(x)),x/Fnorm(x),0))
    if(modality=="All"){Methyl_cell_loading_CH=apply(Methyl_cell_loading_CH,2,function(x)ifelse(is.finite(x/Fnorm(x)),x/Fnorm(x),0))}
    
    }
    
    
    if(modality=="All"|modality=="HiC+CG"){
      
      methyl_res_CG$Methyl_loci_loading=Methyl_loci_loading_CG
      methyl_res_CG$Methyl_cell_loading=Methyl_cell_loading_CG
      if(modality=="All"){ methyl_res_CH$Methyl_loci_loading=Methyl_loci_loading_CH
      methyl_res_CH$Methyl_cell_loading=Methyl_cell_loading_CH}
      
    }

   
    
    
    
    f_old=f_new
    
    
    Theta_list_new=mapply(function(x,y)BTD_combi2(x$A,x$B,x$C,R,y),btd_res,Lr_list,SIMPLIFY = F)
    if(modality=="All"|modality=="HiC+CG"){
      Theta_list_new[['CG']]=methyl_res_CG$Methyl_loci_loading%*%t(methyl_res_CG$Methyl_cell_loading)
      if(modality=="All"){Theta_list_new[['CH']]=methyl_res_CH$Methyl_loci_loading%*%t(methyl_res_CH$Methyl_cell_loading)  }
      
    }
    
    
    
    ls_list=mapply(function(x,y){(sum((x-y)^2))},data,Theta_list_new,SIMPLIFY = F)
    #############
    #update weight
    #############
    
    
    
    f_new=Reduce('+',mapply(function(s,d){1/d*s},ls_list,d_list,SIMPLIFY = F))
    
    
    
    
    
    
    
    
    
    
    
    #####    
    
    
    
    diff=(f_old-f_new)/f_old
    cat(diff)
    
    if((diff)<tol||iter==maxiter){
      if(diff<0){case=1}
      if(diff>=0){case=2}
      break()}
    
  }
  
  print("Rank 1 Muscle over saving results...")
  saveRDS(btd_res,paste0(dir_out,'/HiC_result_iter_',iter,'_rank',k,'.rds'))
  if(modality=="All"|modality=="HiC+CG"){  
    
    saveRDS(methyl_res_CG,paste0(dir_out,'/methyl_CG_result_iter_',iter,'_rank',k,'.rds'))
    if(modality=="All"){ saveRDS(methyl_res_CH,paste0(dir_out,'/methyl_CH_result_iter_',iter,'_rank',k,'.rds')) }}

  
  
  if(case==1){
    rm(Theta_list_new)
    Theta=mapply(function(x,y)BTD_combi2(x$A,x$B,x$C,R,y),btd_res_old,Lr_list,SIMPLIFY = F)
    if(modality=="All"|modality=="HiC+CG"){ 
      Theta[['CG']]=methyl_res_CG_old$Methyl_loci_loading%*%t(methyl_res_CG_old$Methyl_cell_loading)
      if(modality=="All"){ Theta[['CH']]=methyl_res_CH_old$Methyl_loci_loading%*%t(methyl_res_CH_old$Methyl_cell_loading) }
      }
    
    saveRDS(f_old,paste0(dir_out,'/trainerr_',k,'.rds'))
    if(modality=="All"){return(list(btd_res=btd_res_old,methyl_res_CG=methyl_res_CG_old,methyl_res_CH=methyl_res_CH_old,Theta=Theta))}
    if(modality=="HiC"){return(list(btd_res=btd_res_old,Theta=Theta))}
    if(modality=="HiC+CG"){return(list(btd_res=btd_res_old,methyl_res_CG=methyl_res_CG_old,Theta=Theta))}
    
    }
  
  
  if(case==2){saveRDS(f_new,paste0(dir_out,'/trainerr_',k,'.rds'));
    if(modality=="All"){return(list(btd_res=btd_res,methyl_res_CG=methyl_res_CG,methyl_res_CH=methyl_res_CH,Theta=Theta_list_new))}
    if(modality=="HiC"){return(list(btd_res=btd_res,Theta=Theta_list_new))}
    if(modality=="HiC+CG"){return(list(btd_res=btd_res,methyl_res_CG=methyl_res_CG,Theta=Theta_list_new))}
    
    
    }
  
  
} 



