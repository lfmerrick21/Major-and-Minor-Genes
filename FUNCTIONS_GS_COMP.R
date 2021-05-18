calc_missrate <- function(gt_mat)
{
  col_func <- function(gt_col)
  {
    missrate <- sum(is.na(gt_col)) / length(gt_col)
    return(missrate)
  }
  
  missrate_vect <- apply(gt_mat, 2, col_func)
  
  return(missrate_vect)
}

# Calculates the minor allele frequency for every marker in a genotype matrix (coded as c(-1,0,1))
calc_maf_apply <- function(gt_mat, encoding = c(-1, 0, 1))
{
  col_func1 <- function(gt_col)
  {
    allele1_ct <- (sum(gt_col == -1, na.rm = T) * 2) + sum(gt_col == 0, na.rm = T)
    allele2_ct <- (sum(gt_col == 1, na.rm = T) * 2) + sum(gt_col == 0, na.rm = T)
    
    maf <- min(c(allele1_ct, allele2_ct)) / (sum(!is.na(gt_col))*2)
  }
  
  col_func2 <- function(gt_col)
  {
    allele1_ct <- (sum(gt_col == 0, na.rm = T) * 2) + sum(gt_col == 1, na.rm = T)
    allele2_ct <- (sum(gt_col == 2, na.rm = T) * 2) + sum(gt_col == 1, na.rm = T)
    
    maf <- min(c(allele1_ct, allele2_ct)) / (sum(!is.na(gt_col))*2)
  }
  
  if (all(encoding == c(-1, 0, 1)))
  {
    maf_vect <- apply(gt_mat, 2, col_func1)
  } else if (all(encoding == c(0, 1, 2)))
  {
    maf_vect <- apply(gt_mat, 2, col_func2)
  } else{
    print('Encoding not recognized, returning NULL')
    maf_vect <- NULL
  }
  
  return(maf_vect)
}

# This is a function that will split data into a list of k-folds
make_CV_sets <- function(list_length, k = 5)
{
  rand_values <- rnorm(list_length)
  k_quantiles <- quantile(rand_values, 0:k/k)
  k_assign <- cut(rand_values, k_quantiles, include.lowest = T, labels = F)
  
  cv_list <- list()
  for (i in 1:k)
  {
    fold_assignment <- k_assign != i
    cv_list[[i]] <- fold_assignment
  }
  return(cv_list)
}

test_all_models_BLUP_pc_mean_recode <- function(genotypes, phenotype,PCA=NULL,CV=NULL,nIter = 5000, burnIn = 2000, folds = 5)
{
  
  library(BGLR)
  library(rrBLUP)
  library(caret)
  library(tidyr)
  library(dplyr)
  # Make the CV list
  fold_list <- make_CV_sets(length(phenotype[,2]), k = folds)
  
  BGLR_acc_results <- list()
  BGLR_acc_predictions<- list()
  Predictions_ALL<-c()
  for (i in 1:length(fold_list))
  {
    fold_indices <- which(fold_list[[i]])
    
    # Split into training and testing data
    pheno_train <- phenotype[fold_indices,]
    pheno_test=phenotype[-fold_indices,]
    
    # Calculate the GS model using rrBLUP
    myGD_train <- as.matrix(genotypes[fold_indices,])
    myGD_train=apply(myGD_train,2,function(x) recode(x,"0"="-1","1"="0","2"="1"))
    myGD_train=apply(myGD_train,2,as.numeric)
    myGD_test <- as.matrix(genotypes[-fold_indices,])
    myGD_test=apply(myGD_test,2,function(x) recode(x,"0"="-1","1"="0","2"="1"))
    myGD_test=apply(myGD_test,2,as.numeric)
    myY_train <- phenotype[fold_indices,2]
    myY_test <- phenotype[-fold_indices,2]
    myPCA_train <- PCA[fold_indices,]
    myPCA_test <- PCA[-fold_indices,]
    gc()

        if(length(CV)==0){rrBLUP_model <- mixed.solve(y = myY_train,
                                Z = myGD_train)
    pred_effects <- myGD_test %*% rrBLUP_model$u
    fix_effects <- rrBLUP_model$beta
    predictions <- c(pred_effects) + c(fix_effects)
    acc <- cor(predictions, myY_test, use = "pairwise.complete")
    sacc <- cor(predictions, myY_test, use = "pairwise.complete", method = c("spearman"))
    metrics=postResample(pred=predictions,obs=myY_test)
    results=c(ACC=acc,SACC=sacc,metrics)
    prediction=data.frame(phenotype[-fold_indices,],GEBV=predictions,RE=pred_effects,FE=fix_effects)
    
    rrBLUP_model_PC <- mixed.solve(y = myY_train,
                                   Z = myGD_train,
                                   X = myPCA_train)
    pred_effects_PC <- myGD_test %*% rrBLUP_model_PC$u
    fix_effects_PC <- myPCA_test  %*% rrBLUP_model_PC$beta
    predictions_PC <- c(pred_effects_PC) + c(fix_effects_PC)
    acc_PC <- cor(predictions_PC, myY_test, use = "pairwise.complete")
    sacc_PC <- cor(predictions_PC, myY_test, use = "pairwise.complete", method = c("spearman"))
    metrics_PC=postResample(pred=predictions_PC,obs=myY_test)
    results_PC=c(ACC_PC=acc_PC,SACC_PC=sacc_PC,metrics)
    prediction_PC=data.frame(phenotype[-fold_indices,],GEBV=predictions_PC,RE=pred_effects_PC,FE=fix_effects_PC)}else{

    #CV=impute(as.factor(CV))
    if(length(ncol(CV))==0){
      myCV_train <- CV[fold_indices]
      myCV_test <- CV[-fold_indices]
    }else{
      myCV_train <- CV[fold_indices,]
      myCV_test <- CV[-fold_indices,]
    }
    
    fix_train <- as.matrix(myCV_train)
    fix_test  <- as.matrix(myCV_test)
    
    p <- ncol(fix_train)
    XtX <- crossprod(fix_train, fix_train)
    rank.X <- qr(XtX)$rank
    if (rank.X < p) {
      sm2=findLinearCombos(fix_train)$remove
      fix_train=fix_train[,-sm2]
      fix_test=fix_test[,-sm2]}
    
    gc()
    rrBLUP_model <- mixed.solve(y = myY_train,
                                Z = myGD_train,
                                X = fix_train)
    
    pred_effects <- myGD_test %*% rrBLUP_model$u
    fix_effects <- fix_test  %*% rrBLUP_model$beta
    predictions <- c(pred_effects) + c(fix_effects)
    
    acc <- cor(predictions, myY_test, use = "pairwise.complete")
    sacc <- cor(predictions, myY_test, use = "pairwise.complete", method = c("spearman"))
    metrics=postResample(pred=predictions,obs=myY_test)
    results=c(ACC=acc,SACC=sacc,metrics)
    prediction=data.frame(phenotype[-fold_indices,],GEBV=predictions,RE=pred_effects,FE=fix_effects)
    
    fix_train_PC <- as.matrix(cbind(myCV_train,myPCA_train))
    fix_test_PC  <- as.matrix(cbind(myCV_test,myPCA_test))
    
    p <- ncol(fix_train_PC)
    XtX <- crossprod(fix_train_PC, fix_train_PC)
    rank.X <- qr(XtX)$rank
    if (rank.X < p) {
      sm2=findLinearCombos(fix_train_PC)$remove
    fix_train_PC=fix_train_PC[,-sm2]
    fix_test_PC=fix_test_PC[,-sm2]}
    
    rrBLUP_model_PC <- mixed.solve(y = myY_train,
                                   Z = myGD_train,
                                   X = fix_train_PC)
    pred_effects_PC <- myGD_test %*% rrBLUP_model_PC$u
    fix_effects_PC <- fix_test_PC  %*% rrBLUP_model_PC$beta
    predictions_PC <- c(pred_effects_PC) + c(fix_effects_PC)
    acc_PC <- cor(predictions_PC, myY_test, use = "pairwise.complete")
    sacc_PC <- cor(predictions_PC, myY_test, use = "pairwise.complete", method = c("spearman"))
    metrics_PC=postResample(pred=predictions_PC,obs=myY_test)
    results_PC=c(ACC_PC=acc_PC,SACC_PC=sacc_PC,metrics)
    prediction_PC=data.frame(phenotype[-fold_indices,],GEBV=predictions_PC,RE=pred_effects_PC,FE=fix_effects_PC)
    
    }
    Predictions<-cbind(prediction,prediction_PC[,3:5])
    Predictions_ALL=rbind(Predictions_ALL,Predictions)
    BGLR_acc_results[[i]] <- list(results,results_PC)
    #BGLR_acc_predictions[[i]] <- list(results,results_PC,prediction,prediction_PC)
  }
  #, GBLUP_acc, "GBLUP"
  model_vect <- c("Pearson","Spearman","RMSE","R2","MAE", "Pearson_PC","Spearman_PC","RMSE_PC","R2_PC","MAE_PC")
  BGLR_acc_table <- data.frame(matrix(nrow = 0, ncol = 3))
  for (i in 1:length(BGLR_acc_results))
  {
    results_long <- data.frame(rep(i, length(model_vect)), model_vect, unlist(BGLR_acc_results[[i]]))
    BGLR_acc_table <- rbind(BGLR_acc_table, results_long)
  }
  names(BGLR_acc_table) <- c("fold", "model", "r")
  data_wide <- spread(BGLR_acc_table, model, r)
  acc_fold=data.frame(data_wide)
  results=colMeans(acc_fold[2:11], na.rm=TRUE)
  names(Predictions)[6:8]<-c("GEBV_PC","RE_PC","FE_PC")
  results_ALL=list(results,Predictions_ALL)
  return(results_ALL)
}

BLUP_pc_mean_GAGS <- function(genotypes, phenotype,Y=NULL,GM=NULL,GD=NULL,PCA=NULL,CV=NULL,GWAS="BLINK",alpha=0.05, folds = 5)
{
  library(rrBLUP)
  library(tidyr)
  library(Hmisc)
  library(caret)
  library(dplyr)
  library(GAPIT3)
  # Make the CV list
  fold_list <- make_CV_sets(length(phenotype[,2]), k = folds)
  
  BGLR_acc_results <- list()
  BGLR_acc_predictions<- list()
  Predictions_ALL<-list()
  for (i in 1:length(fold_list)){
    fold_indices <- which(fold_list[[i]])
    
    # Split into training and testing data
    # Calculate the GS model using rrBLUP
    myGD_train <- as.matrix(genotypes[fold_indices,])
    myGD_train=apply(myGD_train,2,function(x) recode(x,"0"="-1","1"="0","2"="1"))
    myGD_train=apply(myGD_train,2,as.numeric)
    myGD_test <- as.matrix(genotypes[-fold_indices,])
    myGD_test=apply(myGD_test,2,function(x) recode(x,"0"="-1","1"="0","2"="1"))
    myGD_test=apply(myGD_test,2,as.numeric)
    
    GD_train <- GD[fold_indices,]
    GD_test <- GD[-fold_indices,]
    Y_train <-Y[fold_indices,c(1,2)]
    Y_test <-Y[-fold_indices,c(1,2)]
    myY_train <- phenotype[fold_indices,2]
    myY_test <- phenotype[-fold_indices,2]
    
    GWASR<- GAPIT(Y = Y_train,
                  GD = GD_train,
                  GM = GM,
                  PCA.total=3,
                  model = GWAS,
                  file.output=F)
    
    GWASSM <- which(GWASR$GWAS$P.value <= alpha/length(GWASR$GWAS$P.value))
    if(length(GWASSM)==0){
      acc <- NA
      sacc <- NA
      metrics=c(RMSE=NA,Rsquared=NA,MAE=NA)
      results=c(ACC=acc,SACC=sacc,metrics)
      
      acc_PC <- NA
      sacc_PC <- NA
      metrics_PC=c(RMSE=NA,Rsquared=NA,MAE=NA)
      results_PC=c(ACC_PC=acc_PC,SACC_PC=sacc_PC,metrics)
      
      #prediction=NA
      #prediction_PC=NA
      #Predictions<-c(prediction,prediction_PC)

    }else{
      sm=as.character(GWASR$GWAS[GWASSM,]$SNP)
      
      myCV_train <- myGD_train[,sm]
      myCV_test <- myGD_test[,sm]
      myPCA_train <- PCA[fold_indices,]
      myPCA_test <- PCA[-fold_indices,]
      
      
      fix_train <- as.matrix(myCV_train)
      fix_test  <- as.matrix(myCV_test)
      p <- ncol(fix_train)
      XtX <- crossprod(fix_train, fix_train)
      rank.X <- qr(XtX)$rank
      if (rank.X < p) {
        sm2=findLinearCombos(fix_train)$remove
        fix_train=fix_train[,-sm2]
        fix_test=fix_test[,-sm2]}
      gc()
      rrBLUP_model <- mixed.solve(y = myY_train,
                                  Z = myGD_train,
                                  X = fix_train)
      
      pred_effects <- myGD_test %*% rrBLUP_model$u
      fix_effects <- fix_test  %*% rrBLUP_model$beta
      predictions <- c(pred_effects) + c(fix_effects)
      acc <- cor(predictions, myY_test, use = "pairwise.complete")
      sacc <- cor(predictions, myY_test, use = "pairwise.complete", method = c("spearman"))
      metrics=postResample(pred=predictions,obs=myY_test)
      results=c(ACC=acc,SACC=sacc,metrics)
      #prediction=data.frame(Y_test,GEBV=predictions,RE=pred_effects,FE=fix_effects)
      
      
      
      
      fix_train_PC <- as.matrix(cbind(myCV_train,myPCA_train))
      fix_test_PC  <- as.matrix(cbind(myCV_test,myPCA_test))
      p <- ncol(fix_train_PC)
      XtX <- crossprod(fix_train_PC, fix_train_PC)
      rank.X <- qr(XtX)$rank
      if (rank.X < p) {
        sm2=findLinearCombos(fix_train_PC)$remove
        fix_train_PC=fix_train_PC[,-sm2]
        fix_test_PC=fix_test_PC[,-sm2]}
      gc()
      rrBLUP_model_PC <- mixed.solve(y = myY_train,
                                     Z = myGD_train,
                                     X = fix_train_PC)
      
      pred_effects_PC <- myGD_test %*% rrBLUP_model_PC$u
      fix_effects_PC <- fix_test_PC  %*% rrBLUP_model_PC$beta
      predictions_PC <- c(pred_effects_PC) + c(fix_effects_PC)
      acc_PC <- cor(predictions_PC, myY_test, use = "pairwise.complete")
      sacc_PC <- cor(predictions_PC, myY_test, use = "pairwise.complete", method = c("spearman"))
      metrics_PC=postResample(pred=predictions_PC,obs=myY_test)
      results_PC=c(ACC_PC=acc_PC,SACC_PC=sacc_PC,metrics)
      #prediction_PC=data.frame(Y_test,GEBV=predictions_PC,RE=pred_effects_PC,FE=fix_effects_PC)
     
      #Predictions<-cbind(prediction,prediction_PC[,3:5])
      
    }
    
    #Predictions_ALL[[i]]=Predictions
    BGLR_acc_results[[i]] <- list(results,results_PC)
    #BGLR_acc_predictions[[i]] <- list(results,results_PC,prediction,prediction_PC)
  }
  #, GBLUP_acc, "GBLUP"
  model_vect <- c("Pearson","Spearman","RMSE","R2","MAE", "Pearson_PC","Spearman_PC","RMSE_PC","R2_PC","MAE_PC")
  BGLR_acc_table <- data.frame(matrix(nrow = 0, ncol = 3))
  for (i in 1:length(BGLR_acc_results))
  {
    results_long <- data.frame(rep(i, length(model_vect)), model_vect, unlist(BGLR_acc_results[[i]]))
    BGLR_acc_table <- rbind(BGLR_acc_table, results_long)
  }
  names(BGLR_acc_table) <- c("fold", "model", "r")
  data_wide <- spread(BGLR_acc_table, model, r)
  acc_fold=data.frame(data_wide)
  results=colMeans(acc_fold[2:11], na.rm=TRUE)
  #names(Predictions)[6:8]<-c("GEBV_PC","RE_PC","FE_PC")
  #results_ALL=list(results,Predictions_ALL)
  #return(results_ALL)
  return(results)
}

BLUP_pc_mean_GAGS_10 <- function(genotypes, phenotype,Y=NULL,GM=NULL,GD=NULL,PCA=NULL,CV=NULL,GWAS="BLINK",alpha=0.05, folds = 5)
{
  library(rrBLUP)
  library(tidyr)
  library(Hmisc)
  library(caret)
  # Make the CV list
  fold_list <- make_CV_sets(length(phenotype[,2]), k = folds)
  
  BGLR_acc_results <- list()
  BGLR_acc_predictions<- list()
  Predictions_ALL<-list()
  for (i in 1:length(fold_list)){
    fold_indices <- which(fold_list[[i]])
    
    # Split into training and testing data
    # Calculate the GS model using rrBLUP
    myGD_train <- as.matrix(genotypes[fold_indices,])
    myGD_train=apply(myGD_train,2,function(x) recode(x,"0"="-1","1"="0","2"="1"))
    myGD_train=apply(myGD_train,2,as.numeric)
    myGD_test <- as.matrix(genotypes[-fold_indices,])
    myGD_test=apply(myGD_test,2,function(x) recode(x,"0"="-1","1"="0","2"="1"))
    myGD_test=apply(myGD_test,2,as.numeric)
    
    GD_train <- GD[fold_indices,]
    GD_test <- GD[-fold_indices,]
    Y_train <-Y[fold_indices,c(1,2)]
    Y_test <-Y[-fold_indices,c(1,2)]
    myY_train <- phenotype[fold_indices,2]
    myY_test <- phenotype[-fold_indices,2]
    
    GWASR<- GAPIT(Y = Y_train,
                  GD = GD_train,
                  GM = GM,
                  PCA.total=3,
                  model = GWAS,
                  file.output=F)
    
    top10=GWASR$GWAS[order(GWASR$GWAS$P.value),]
    GWASSM=top10[1:10,]$SNP
    myCV_train <- myGD_train[,GWASSM]
    myCV_test <- myGD_test[,GWASSM]
    myPCA_train <- PCA[fold_indices,]
    myPCA_test <- PCA[-fold_indices,]
    
    fix_train <- as.matrix(myCV_train)
    fix_test  <- as.matrix(myCV_test)
    p <- ncol(fix_train)
    XtX <- crossprod(fix_train, fix_train)
    rank.X <- qr(XtX)$rank
    if (rank.X < p) {
      sm2=findLinearCombos(fix_train)$remove
      fix_train=fix_train[,-sm2]
      fix_test=fix_test[,-sm2]}
    gc()
    rrBLUP_model <- mixed.solve(y = myY_train,
                                Z = myGD_train,
                                X = fix_train)
    
    pred_effects <- myGD_test %*% rrBLUP_model$u
    fix_effects <- fix_test  %*% rrBLUP_model$beta
    predictions <- c(pred_effects) + c(fix_effects)
    acc <- cor(predictions, myY_test, use = "pairwise.complete")
    sacc <- cor(predictions, myY_test, use = "pairwise.complete", method = c("spearman"))
    metrics=postResample(pred=predictions,obs=myY_test)
    results=c(ACC=acc,SACC=sacc,metrics)
    prediction=data.frame(Y_test,GEBV=predictions,RE=pred_effects,FE=fix_effects)
    
    
    
    
    fix_train_PC <- as.matrix(cbind(myCV_train,myPCA_train))
    fix_test_PC  <- as.matrix(cbind(myCV_test,myPCA_test))
    p <- ncol(fix_train_PC)
    XtX <- crossprod(fix_train_PC, fix_train_PC)
    rank.X <- qr(XtX)$rank
    if (rank.X < p) {
      sm2=findLinearCombos(fix_train_PC)$remove
      fix_train_PC=fix_train_PC[,-sm2]
      fix_test_PC=fix_test_PC[,-sm2]}
    gc()
    rrBLUP_model_PC <- mixed.solve(y = myY_train,
                                   Z = myGD_train,
                                   X = fix_train_PC)
    
    pred_effects_PC <- myGD_test %*% rrBLUP_model_PC$u
    fix_effects_PC <- fix_test_PC  %*% rrBLUP_model_PC$beta
    predictions_PC <- c(pred_effects_PC) + c(fix_effects_PC)
    acc_PC <- cor(predictions_PC, myY_test, use = "pairwise.complete")
    sacc_PC <- cor(predictions_PC, myY_test, use = "pairwise.complete", method = c("spearman"))
    metrics_PC=postResample(pred=predictions_PC,obs=myY_test)
    results_PC=c(ACC_PC=acc_PC,SACC_PC=sacc_PC,metrics)
    prediction_PC=data.frame(Y_test,GEBV=predictions_PC,RE=pred_effects_PC,FE=fix_effects_PC)
    
    
    
  
  Predictions<-cbind(prediction,prediction_PC[,3:5])
  Predictions_ALL[[i]]=Predictions
  BGLR_acc_results[[i]] <- list(results,results_PC)
  #BGLR_acc_predictions[[i]] <- list(results,results_PC,prediction,prediction_PC)
  }
  #, GBLUP_acc, "GBLUP"
  model_vect <- c("Pearson","Spearman","RMSE","R2","MAE", "Pearson_PC","Spearman_PC","RMSE_PC","R2_PC","MAE_PC")
  BGLR_acc_table <- data.frame(matrix(nrow = 0, ncol = 3))
  for (i in 1:length(BGLR_acc_results))
  {
    results_long <- data.frame(rep(i, length(model_vect)), model_vect, unlist(BGLR_acc_results[[i]]))
    BGLR_acc_table <- rbind(BGLR_acc_table, results_long)
  }
  names(BGLR_acc_table) <- c("fold", "model", "r")
  data_wide <- spread(BGLR_acc_table, model, r)
  acc_fold=data.frame(data_wide)
  results=colMeans(acc_fold[2:11], na.rm=TRUE)
  names(Predictions)[6:8]<-c("GEBV_PC","RE_PC","FE_PC")
  results_ALL=list(results,Predictions_ALL)
  return(results_ALL)
}

MAS_GWAS_CV_Mean <- function(genotypes, phenotype,Y=NULL,GM=NULL,GD=NULL,PCA=NULL,CV=NULL,GWAS="BLINK",alpha=0.05, folds = 5)
{
  library(rrBLUP)
  library(tidyr)
  library(Hmisc)
  library(caret)
  library(dplyr)
  # Make the CV list
  fold_list <- make_CV_sets(length(phenotype[,2]), k = folds)
  
  BGLR_acc_results <- list()
  BGLR_acc_predictions<- list()
  Predictions_ALL<-list()
  for (i in 1:length(fold_list)){
    fold_indices <- which(fold_list[[i]])
    
    # Split into training and testing data
    # Calculate the GS model using rrBLUP
    myGD_train <- as.matrix(genotypes[fold_indices,])
    myGD_train=apply(myGD_train,2,function(x) recode(x,"0"="-1","1"="0","2"="1"))
    myGD_train=apply(myGD_train,2,as.numeric)
    myGD_test <- as.matrix(genotypes[-fold_indices,])
    myGD_test=apply(myGD_test,2,function(x) recode(x,"0"="-1","1"="0","2"="1"))
    myGD_test=apply(myGD_test,2,as.numeric)
    
    GD_train <- GD[fold_indices,]
    GD_test <- GD[-fold_indices,]
    Y_train <-Y[fold_indices,c(1,2)]
    Y_test <-Y[-fold_indices,c(1,2)]
    
    pheno_train <- phenotype[,2]
    pheno_train[-fold_indices] <- NA
    
    myY_train <- phenotype[fold_indices,2]
    myY_test <- phenotype[-fold_indices,2]
    
    GWASR<- GAPIT(Y = Y_train,
                  GD = GD_train,
                  GM = GM,
                  PCA.total=3,
                  model = GWAS,
                  file.output=F)
    
    GWASSM <- which(GWASR$GWAS$P.value <= alpha/length(GWASR$GWAS$P.value))
    if(length(GWASSM)==0){
      acc <- NA
      sacc <- NA
      metrics=c(RMSE=NA,Rsquared=NA,MAE=NA)
      results=c(ACC=acc,SACC=sacc,metrics)
      
      acc_PC <- NA
      sacc_PC <- NA
      metrics_PC=c(RMSE=NA,Rsquared=NA,MAE=NA)
      results_PC=c(ACC_PC=acc_PC,SACC_PC=sacc_PC,metrics)
      
      prediction=NA
      prediction_PC=NA
      Predictions<-c(prediction,prediction_PC)
    }else{
      sm=as.character(GWASR$GWAS[GWASSM,]$SNP)
      
      CV <- genotypes[,sm]
      myPCA_train <- PCA
      myPCA_test <- PCA
      
      
      MAS_train <- data.frame(CV)
      MAS_test  <- data.frame(CV)
      
      GLM_data <- data.frame(pheno_train, MAS_train)

      names(GLM_data)[1] <- "Y"
      #Linear model to calculate effects
      #You can run all signficant markers at once to see cumulative effect
      MAS_model <- lm(Y ~ ., data = GLM_data)
      MAS_pred <- predict(MAS_model, MAS_test)
      
      acc <- cor(MAS_pred[-fold_indices], myY_test, use = "pairwise.complete")
      sacc <- cor(MAS_pred[-fold_indices], myY_test, use = "pairwise.complete", method = c("spearman"))
      metrics=postResample(pred=MAS_pred[-fold_indices],obs=myY_test)
      results=c(ACC=acc,SACC=sacc,metrics)
      #prediction=data.frame(Y_test,GEBV=MAS_pred[-fold_indices])
      
      
      MAS_train_PC <- data.frame(cbind(MAS_train,myPCA_train))
      MAS_test_PC  <- data.frame(cbind(MAS_test,myPCA_test))
      
      GLM_data_PC <- data.frame(pheno_train, MAS_train_PC)
 
      names(GLM_data_PC)[1] <- "Y"
      #Linear model to calculate effects
      #You can run all signficant markers at once to see cumulative effect
      MAS_model_PC <- lm(Y ~ ., data = GLM_data_PC)
      MAS_pred_PC <- predict(MAS_model_PC, MAS_test_PC)
      
      acc_PC <- cor(MAS_pred_PC[-fold_indices], myY_test, use = "pairwise.complete")
      sacc_PC <- cor(MAS_pred_PC[-fold_indices], myY_test, use = "pairwise.complete", method = c("spearman"))
      metrics_PC=postResample(pred=MAS_pred_PC[-fold_indices],obs=myY_test)
      results_PC=c(ACC_PC=acc_PC,SACC_PC=sacc_PC,metrics)
      #prediction_PC=data.frame(Y_test,GEBV=MAS_pred_PC[-fold_indices])
      #Predictions<-cbind(prediction,prediction_PC[,3])
      
    }
    
    #Predictions_ALL[[i]]=Predictions
    BGLR_acc_results[[i]] <- list(results,results_PC)
    #BGLR_acc_predictions[[i]] <- list(results,results_PC,prediction,prediction_PC)
  }
  #, GBLUP_acc, "GBLUP"
  model_vect <- c("Pearson","Spearman","RMSE","R2","MAE", "Pearson_PC","Spearman_PC","RMSE_PC","R2_PC","MAE_PC")
  BGLR_acc_table <- data.frame(matrix(nrow = 0, ncol = 3))
  for (i in 1:length(BGLR_acc_results))
  {
    results_long <- data.frame(rep(i, length(model_vect)), model_vect, unlist(BGLR_acc_results[[i]]))
    BGLR_acc_table <- rbind(BGLR_acc_table, results_long)
  }
  names(BGLR_acc_table) <- c("fold", "model", "r")
  data_wide <- spread(BGLR_acc_table, model, r)
  acc_fold=data.frame(data_wide)
  results=colMeans(acc_fold[2:11], na.rm=TRUE)
  #names(Predictions)[4]<-c("GEBV_PC")
  #results_ALL=list(results,Predictions_ALL)
  #return(results_ALL)
  return(results)
}

MAS_GWAS_10_CV_Mean <- function(genotypes, phenotype,Y=NULL,GM=NULL,GD=NULL,PCA=NULL,CV=NULL,GWAS="BLINK",alpha=0.05, folds = 5)
{
  library(rrBLUP)
  library(tidyr)
  library(Hmisc)
  library(caret)
  library(dplyr)
  # Make the CV list
  fold_list <- make_CV_sets(length(phenotype[,2]), k = folds)
  
  BGLR_acc_results <- list()
  BGLR_acc_predictions<- list()
  Predictions_ALL<-list()
  for (i in 1:length(fold_list)){
    fold_indices <- which(fold_list[[i]])
    
    # Split into training and testing data
    # Calculate the GS model using rrBLUP
    myGD_train <- as.matrix(genotypes[fold_indices,])
    myGD_train=apply(myGD_train,2,function(x) recode(x,"0"="-1","1"="0","2"="1"))
    myGD_train=apply(myGD_train,2,as.numeric)
    myGD_test <- as.matrix(genotypes[-fold_indices,])
    myGD_test=apply(myGD_test,2,function(x) recode(x,"0"="-1","1"="0","2"="1"))
    myGD_test=apply(myGD_test,2,as.numeric)
    
    GD_train <- GD[fold_indices,]
    GD_test <- GD[-fold_indices,]
    Y_train <-Y[fold_indices,c(1,2)]
    Y_test <-Y[-fold_indices,c(1,2)]
    
    pheno_train <- phenotype[,2]
    pheno_train[-fold_indices] <- NA
    
    myY_train <- phenotype[fold_indices,2]
    myY_test <- phenotype[-fold_indices,2]
    
    GWASR<- GAPIT(Y = Y_train,
                  GD = GD_train,
                  GM = GM,
                  PCA.total=3,
                  model = GWAS,
                  file.output=F)
    
    top10=GWASR$GWAS[order(GWASR$GWAS$P.value),]
    GWASSM=top10[1:10,]$SNP
    
    
    CV <- genotypes[,GWASSM]
    myPCA_train <- PCA
    myPCA_test <- PCA
    
    
    MAS_train <- data.frame(CV)
    MAS_test  <- data.frame(CV)
    
    GLM_data <- data.frame(pheno_train, MAS_train)

    names(GLM_data)[1] <- "Y"
    #Linear model to calculate effects
    #You can run all signficant markers at once to see cumulative effect
    MAS_model <- lm(Y ~ ., data = GLM_data)
    MAS_pred <- predict(MAS_model, MAS_test)
    
    acc <- cor(MAS_pred[-fold_indices], myY_test, use = "pairwise.complete")
    sacc <- cor(MAS_pred[-fold_indices], myY_test, use = "pairwise.complete", method = c("spearman"))
    metrics=postResample(pred=MAS_pred[-fold_indices],obs=myY_test)
    results=c(ACC=acc,SACC=sacc,metrics)
    prediction=data.frame(Y_test,GEBV=MAS_pred[-fold_indices])
    
    
    
    
    MAS_train_PC <- data.frame(cbind(MAS_train,myPCA_train))
    MAS_test_PC  <- data.frame(cbind(MAS_test,myPCA_test))
    
    GLM_data_PC <- data.frame(pheno_train, MAS_train_PC)

    names(GLM_data_PC)[1] <- "Y"
    #Linear model to calculate effects
    #You can run all signficant markers at once to see cumulative effect
    MAS_model_PC <- lm(Y ~ ., data = GLM_data_PC)
    MAS_pred_PC <- predict(MAS_model_PC, MAS_test_PC)
    
    acc_PC <- cor(MAS_pred_PC[-fold_indices], myY_test, use = "pairwise.complete")
    sacc_PC <- cor(MAS_pred_PC[-fold_indices], myY_test, use = "pairwise.complete", method = c("spearman"))
    metrics_PC=postResample(pred=MAS_pred_PC[-fold_indices],obs=myY_test)
    results_PC=c(ACC_PC=acc_PC,SACC_PC=sacc_PC,metrics)
    prediction_PC=data.frame(Y_test,GEBV=MAS_pred_PC[-fold_indices])
    
    
    Predictions<-cbind(prediction,prediction_PC[,3])
    Predictions_ALL[[i]]=Predictions
    BGLR_acc_results[[i]] <- list(results,results_PC)
    #BGLR_acc_predictions[[i]] <- list(results,results_PC,prediction,prediction_PC)
  }
  #, GBLUP_acc, "GBLUP"
  model_vect <- c("Pearson","Spearman","RMSE","R2","MAE", "Pearson_PC","Spearman_PC","RMSE_PC","R2_PC","MAE_PC")
  BGLR_acc_table <- data.frame(matrix(nrow = 0, ncol = 3))
  for (i in 1:length(BGLR_acc_results))
  {
    results_long <- data.frame(rep(i, length(model_vect)), model_vect, unlist(BGLR_acc_results[[i]]))
    BGLR_acc_table <- rbind(BGLR_acc_table, results_long)
  }
  names(BGLR_acc_table) <- c("fold", "model", "r")
  data_wide <- spread(BGLR_acc_table, model, r)
  acc_fold=data.frame(data_wide)
  results=colMeans(acc_fold[2:11], na.rm=TRUE)
  names(Predictions)[4]<-c("GEBV_PC")
  results_ALL=list(results,Predictions_ALL)
  return(results_ALL)
}


MAS_CV_Mean <- function(phenotype,PCA=NULL,CV=NULL,nIter = 5000, burnIn = 2000, folds = 5)
{
  
  library(BGLR)
  library(rrBLUP)
  library(caret)
  library(tidyr)
  library(dplyr)
  # Make the CV list
  fold_list <- make_CV_sets(length(phenotype[,2]), k = folds)
  
  BGLR_acc_results <- list()
  BGLR_acc_predictions<- list()
  Predictions_ALL<-list()
  for (i in 1:length(fold_list))
  {
    fold_indices <- which(fold_list[[i]])
    
    pheno_train <- phenotype[,2]
    pheno_train[-fold_indices] <- NA
    # Calculate the GS model using rrBLUP

    myY_train <- phenotype[fold_indices,2]
    myY_test <- phenotype[-fold_indices,2]

    
    # Calculate the GS model using rrBLUP
    myY_train <- phenotype[fold_indices,2]
    myY_test <- phenotype[-fold_indices,2]
    myPCA_train <- PCA
    myPCA_test <- PCA
    
    gc()

    #CV=impute(as.factor(CV))
    #if(length(ncol(CV))==0){
      #myCV_train <- CV[fold_indices]
      #myCV_test <- CV[-fold_indices]
    #}else{
      #myCV_train <- CV[fold_indices,]
      #myCV_test <- CV[-fold_indices,]
    #}
    MAS_train <- data.frame(CV)
    MAS_test  <- data.frame(CV)
    
    GLM_data <- data.frame(pheno_train, MAS_train)

    names(GLM_data)[1] <- "Y"
    #Linear model to calculate effects
    #You can run all signficant markers at once to see cumulative effect
    MAS_model <- lm(Y ~ ., data = GLM_data)
    MAS_pred <- predict(MAS_model,MAS_test)
    length(MAS_pred)
    length(phenotype[-fold_indices,2])
    length(MAS_pred[-fold_indices])

    acc <- cor(MAS_pred[-fold_indices], myY_test, use = "pairwise.complete")
    sacc <- cor(MAS_pred[-fold_indices], myY_test, use = "pairwise.complete", method = c("spearman"))
    metrics=postResample(pred=MAS_pred[-fold_indices],obs=myY_test)
    results=c(ACC=acc,SACC=sacc,metrics)
    prediction=data.frame(phenotype[-fold_indices,],GEBV=MAS_pred[-fold_indices])
  
    
    MAS_train_PC <- data.frame(cbind(MAS_train,myPCA_train))
    MAS_test_PC  <- data.frame(cbind(MAS_test,myPCA_test))
    
    GLM_data_PC <- data.frame(pheno_train, MAS_train_PC)

    names(GLM_data_PC)[1] <- "Y"
    #Linear model to calculate effects
    #You can run all signficant markers at once to see cumulative effect
    MAS_model_PC <- lm(Y ~ ., data = GLM_data_PC)
    MAS_pred_PC <- predict(MAS_model_PC, MAS_test_PC)
    
    acc_PC <- cor(MAS_pred_PC[-fold_indices], myY_test, use = "pairwise.complete")
    sacc_PC <- cor(MAS_pred_PC[-fold_indices], myY_test, use = "pairwise.complete", method = c("spearman"))
    metrics_PC=postResample(pred=MAS_pred_PC[-fold_indices],obs=myY_test)
    results_PC=c(ACC_PC=acc_PC,SACC_PC=sacc_PC,metrics)
    prediction_PC=data.frame(phenotype[-fold_indices,],GEBV=MAS_pred_PC[-fold_indices])
    
    
    Predictions<-cbind(prediction,prediction_PC[,3])
    Predictions_ALL[[i]]=Predictions
    BGLR_acc_results[[i]] <- list(results,results_PC)
    #BGLR_acc_predictions[[i]] <- list(results,results_PC,prediction,prediction_PC)
  }
  #, GBLUP_acc, "GBLUP"
  model_vect <- c("Pearson","Spearman","RMSE","R2","MAE", "Pearson_PC","Spearman_PC","RMSE_PC","R2_PC","MAE_PC")
  BGLR_acc_table <- data.frame(matrix(nrow = 0, ncol = 3))
  for (i in 1:length(BGLR_acc_results))
  {
    results_long <- data.frame(rep(i, length(model_vect)), model_vect, unlist(BGLR_acc_results[[i]]))
    BGLR_acc_table <- rbind(BGLR_acc_table, results_long)
  }
  names(BGLR_acc_table) <- c("fold", "model", "r")
  data_wide <- spread(BGLR_acc_table, model, r)
  acc_fold=data.frame(data_wide)
  results=colMeans(acc_fold[2:11], na.rm=TRUE)
  names(Predictions)[4]<-c("GEBV_PC")
  results_ALL=list(results,Predictions_ALL)
  return(results_ALL)
}



test_all_models_BLUP_vs_pc_mean_recode <- function(train_genotypes, train_phenotype,train_PCA=NULL,train_CV=NULL,test_genotypes, test_phenotype,test_PCA=NULL,test_CV=NULL)
{
  library(BGLR)
  library(tidyr)
  library(rrBLUP)
  library(caret)
  library(dplyr)
  
  # Calculate the GS model using rrBLUP
  myGD_train <- as.matrix(train_genotypes)
  myGD_train=apply(myGD_train,2,function(x) recode(x,"0"="-1","1"="0","2"="1"))
  myGD_train=apply(myGD_train,2,as.numeric)
  
  myGD_test <- as.matrix(test_genotypes)
  myGD_test=apply(myGD_test,2,function(x) recode(x,"0"="-1","1"="0","2"="1"))
  myGD_test=apply(myGD_test,2,as.numeric)
  
  myY_train <- train_phenotype[,2]
  myY_test <- test_phenotype[,2]
  myPCA_train <- train_PCA
  myPCA_test <- test_PCA
  
  gc()
  if(length(train_CV)==0){rrBLUP_model <- mixed.solve(y = myY_train,
                                                      Z = myGD_train)
  pred_effects <- myGD_test %*% rrBLUP_model$u
  fix_effects <- rrBLUP_model$beta
  predictions <- c(pred_effects) + c(fix_effects)
  acc <- cor(predictions, myY_test, use = "pairwise.complete")
  sacc <- cor(predictions, myY_test, use = "pairwise.complete", method = c("spearman"))
  metrics=postResample(pred=predictions,obs=myY_test)
  results=c(ACC=acc,SACC=sacc,metrics)
  prediction=data.frame(test_phenotype,GEBV=predictions,RE=pred_effects,FE=fix_effects)
  
  rrBLUP_model_PC <- mixed.solve(y = myY_train,
                                 Z = myGD_train,
                                 X = myPCA_train)
  pred_effects_PC <- myGD_test %*% rrBLUP_model_PC$u
  fix_effects_PC <- myPCA_test  %*% rrBLUP_model_PC$beta
  predictions_PC <- c(pred_effects_PC) + c(fix_effects_PC)
  acc_PC <- cor(predictions_PC, myY_test, use = "pairwise.complete")
  sacc_PC <- cor(predictions_PC, myY_test, use = "pairwise.complete", method = c("spearman"))
  metrics_PC=postResample(pred=predictions_PC,obs=myY_test)
  results_PC=c(ACC_PC=acc_PC,SACC_PC=sacc_PC,metrics)
  prediction_PC=data.frame(test_phenotype,GEBV=predictions_PC,RE=pred_effects_PC,FE=fix_effects_PC)
  }else{
    fix_train <- as.matrix(train_CV)
    fix_test  <- as.matrix(test_CV)
    
    p <- ncol(fix_train)
    XtX <- crossprod(fix_train, fix_train)
    rank.X <- qr(XtX)$rank
    if (rank.X < p) {
      sm2=findLinearCombos(fix_train)$remove
      fix_train=fix_train[,-sm2]
      fix_test=fix_test[,-sm2]}
    
    gc()
    rrBLUP_model <- mixed.solve(y = myY_train,
                                Z = myGD_train,
                                X = fix_train)
    
    pred_effects <- myGD_test %*% rrBLUP_model$u
    fix_effects <- fix_test  %*% rrBLUP_model$beta
    predictions <- c(pred_effects) + c(fix_effects)
    
    acc <- cor(predictions, myY_test, use = "pairwise.complete")
    sacc <- cor(predictions, myY_test, use = "pairwise.complete", method = c("spearman"))
    metrics=postResample(pred=predictions,obs=myY_test)
    results=c(ACC=acc,SACC=sacc,metrics)
    prediction=data.frame(test_phenotype,GEBV=predictions,RE=pred_effects,FE=fix_effects)
    
    fix_train_PC <- as.matrix(cbind(train_CV,myPCA_train))
    fix_test_PC  <- as.matrix(cbind(test_CV,myPCA_test))
    
    p <- ncol(fix_train_PC)
    XtX <- crossprod(fix_train_PC, fix_train_PC)
    rank.X <- qr(XtX)$rank
    if (rank.X < p) {
      sm2=findLinearCombos(fix_train_PC)$remove
      fix_train_PC=fix_train_PC[,-sm2]
      fix_test_PC=fix_test_PC[,-sm2]}
    
    rrBLUP_model_PC <- mixed.solve(y = myY_train,
                                   Z = myGD_train,
                                   X = fix_train_PC)
    pred_effects_PC <- myGD_test %*% rrBLUP_model_PC$u
    fix_effects_PC <- fix_test_PC  %*% rrBLUP_model_PC$beta
    predictions_PC <- c(pred_effects_PC) + c(fix_effects_PC)
    acc_PC <- cor(predictions_PC, myY_test, use = "pairwise.complete")
    sacc_PC <- cor(predictions_PC, myY_test, use = "pairwise.complete", method = c("spearman"))
    metrics_PC=postResample(pred=predictions_PC,obs=myY_test)
    results_PC=c(ACC_PC=acc_PC,SACC_PC=sacc_PC,metrics)
    prediction_PC=data.frame(test_phenotype,GEBV=predictions_PC,RE=pred_effects_PC,FE=fix_effects_PC)
    
    
  }
  
  Accuracy=c(results,results_PC)
  names(Accuracy)<-c("Pearson","Spearman","RMSE","R2","MAE", "Pearson_PC","Spearman_PC","RMSE_PC","R2_PC","MAE_PC")
  Predictions<-cbind(prediction,prediction_PC[,3:5])
  names(Predictions)[6:8]<-c("GEBV_PC","RE_PC","FE_PC")
  results_ALL=list(Accuracy,Predictions)
  return(results_ALL)
}

BLUP_vs_pc_mean_GAGS <- function(train_genotypes, train_phenotype,train_Y=NULL,train_GM=NULL,train_GD=NULL,train_PCA=NULL,test_genotypes,test_phenotype,test_PCA=NULL,GWAS="BLINK",alpha=0.05, folds = 5)
{
  library(rrBLUP)
  library(tidyr)
  library(Hmisc)
  library(caret)
  
  myGD_train <- as.matrix(train_genotypes)
  myGD_train=apply(myGD_train,2,function(x) recode(x,"0"="-1","1"="0","2"="1"))
  myGD_train=apply(myGD_train,2,as.numeric)
  myGD_test <- as.matrix(test_genotypes)
  myGD_test=apply(myGD_test,2,function(x) recode(x,"0"="-1","1"="0","2"="1"))
  myGD_test=apply(myGD_test,2,as.numeric)
  
  myY_train <- train_phenotype[,2]
  myY_test <- test_phenotype[,2]
  
  GWASR<- GAPIT(Y = train_Y,
                GD = train_GD,
                GM = train_GM,
                PCA.total=3,
                model = GWAS,
                file.output=F)
  
  GWASSM <- which(GWASR$GWAS$P.value <= alpha/length(GWASR$GWAS$P.value))
  if(length(GWASSM)==0){
    acc <- NA
    sacc <- NA
    metrics=c(RMSE=NA,Rsquared=NA,MAE=NA)
    results=c(ACC=acc,SACC=sacc,metrics)
    prediction=NA
    
    acc_PC <- NA
    sacc_PC <- NA
    metrics_PC=c(RMSE=NA,Rsquared=NA,MAE=NA)
    results_PC=c(ACC_PC=acc_PC,SACC_PC=sacc_PC,metrics)
    
    prediction_PC=NA
    
    Accuracy=c(results,results_PC)
    
    names(Accuracy)<-c("Pearson","Spearman","RMSE","R2","MAE", "Pearson_PC","Spearman_PC","RMSE_PC","R2_PC","MAE_PC")
    Predictions<-cbind(prediction,prediction_PC)
    
    results_ALL=list(Accuracy,Predictions)
  }else{
    sm=as.character(GWASR$GWAS[GWASSM,]$SNP)
    
    myCV_train <- myGD_train[,sm]
    myCV_test <- myGD_test[,sm]
    myPCA_train <- train_PCA
    myPCA_test <- test_PCA
    
    
    fix_train <- as.matrix(myCV_train)
    fix_test  <- as.matrix(myCV_test)
    p <- ncol(fix_train)
    XtX <- crossprod(fix_train, fix_train)
    rank.X <- qr(XtX)$rank
    if (rank.X < p) {
      sm2=findLinearCombos(fix_train)$remove
      fix_train=fix_train[,-sm2]
      fix_test=fix_test[,-sm2]}
    gc()
    rrBLUP_model <- mixed.solve(y = myY_train,
                                Z = myGD_train,
                                X = fix_train)
    
    pred_effects <- myGD_test %*% rrBLUP_model$u
    fix_effects <- fix_test  %*% rrBLUP_model$beta
    predictions <- c(pred_effects) + c(fix_effects)
    acc <- cor(predictions, myY_test, use = "pairwise.complete")
    sacc <- cor(predictions, myY_test, use = "pairwise.complete", method = c("spearman"))
    metrics=postResample(pred=predictions,obs=myY_test)
    results=c(ACC=acc,SACC=sacc,metrics)
    prediction=data.frame(test_phenotype,GEBV=predictions,RE=pred_effects,FE=fix_effects)
    
    
    
    
    fix_train_PC <- as.matrix(cbind(myCV_train,myPCA_train))
    fix_test_PC  <- as.matrix(cbind(myCV_test,myPCA_test))
    p <- ncol(fix_train_PC)
    XtX <- crossprod(fix_train_PC, fix_train_PC)
    rank.X <- qr(XtX)$rank
    if (rank.X < p) {
      sm2=findLinearCombos(fix_train_PC)$remove
      fix_train_PC=fix_train_PC[,-sm2]
      fix_test_PC=fix_test_PC[,-sm2]}
    gc()
    rrBLUP_model_PC <- mixed.solve(y = myY_train,
                                   Z = myGD_train,
                                   X = fix_train_PC)
    
    pred_effects_PC <- myGD_test %*% rrBLUP_model_PC$u
    fix_effects_PC <- fix_test_PC  %*% rrBLUP_model_PC$beta
    predictions_PC <- c(pred_effects_PC) + c(fix_effects_PC)
    acc_PC <- cor(predictions_PC, myY_test, use = "pairwise.complete")
    sacc_PC <- cor(predictions_PC, myY_test, use = "pairwise.complete", method = c("spearman"))
    metrics_PC=postResample(pred=predictions_PC,obs=myY_test)
    results_PC=c(ACC_PC=acc_PC,SACC_PC=sacc_PC,metrics)
    prediction_PC=data.frame(test_phenotype,GEBV=predictions_PC,RE=pred_effects_PC,FE=fix_effects_PC)
    Accuracy=c(results,results_PC)
    names(Accuracy)<-c("Pearson","Spearman","RMSE","R2","MAE", "Pearson_PC","Spearman_PC","RMSE_PC","R2_PC","MAE_PC")
    Predictions<-cbind(prediction,prediction_PC[,3:5])
    names(Predictions)[6:8]<-c("GEBV_PC","RE_PC","FE_PC")
    results_ALL=list(Accuracy,Predictions)
  }
  
  return(results_ALL)
}

BLUP_vs_pc_mean_GAGS_10 <- function(train_genotypes, train_phenotype,train_Y=NULL,train_GM=NULL,train_GD=NULL,train_PCA=NULL,test_genotypes,test_phenotype,test_PCA=NULL,GWAS="BLINK",alpha=0.05, folds = 5)
{
  library(rrBLUP)
  library(tidyr)
  library(Hmisc)
  library(caret)
  
  myGD_train <- as.matrix(train_genotypes)
  myGD_train=apply(myGD_train,2,function(x) recode(x,"0"="-1","1"="0","2"="1"))
  myGD_train=apply(myGD_train,2,as.numeric)
  myGD_test <- as.matrix(test_genotypes)
  myGD_test=apply(myGD_test,2,function(x) recode(x,"0"="-1","1"="0","2"="1"))
  myGD_test=apply(myGD_test,2,as.numeric)
  
  myY_train <- train_phenotype[,2]
  myY_test <- test_phenotype[,2]
  
  GWASR<- GAPIT(Y = train_Y,
                GD = train_GD,
                GM = train_GM,
                PCA.total=3,
                model = GWAS,
                file.output=F)
  
  top10=GWASR$GWAS[order(GWASR$GWAS$P.value),]
  GWASSM=top10[1:10,]$SNP

  
  myCV_train <- myGD_train[,GWASSM]
  myCV_test <- myGD_test[,GWASSM]
  myPCA_train <- train_PCA
  myPCA_test <- test_PCA
  
  
  fix_train <- as.matrix(myCV_train)
  fix_test  <- as.matrix(myCV_test)
  p <- ncol(fix_train)
  XtX <- crossprod(fix_train, fix_train)
  rank.X <- qr(XtX)$rank
  if (rank.X < p) {
    sm2=findLinearCombos(fix_train)$remove
    fix_train=fix_train[,-sm2]
    fix_test=fix_test[,-sm2]}
  gc()
  rrBLUP_model <- mixed.solve(y = myY_train,
                              Z = myGD_train,
                              X = fix_train)
  
  pred_effects <- myGD_test %*% rrBLUP_model$u
  fix_effects <- fix_test  %*% rrBLUP_model$beta
  predictions <- c(pred_effects) + c(fix_effects)
  acc <- cor(predictions, myY_test, use = "pairwise.complete")
  sacc <- cor(predictions, myY_test, use = "pairwise.complete", method = c("spearman"))
  metrics=postResample(pred=predictions,obs=myY_test)
  results=c(ACC=acc,SACC=sacc,metrics)
  prediction=data.frame(test_phenotype,GEBV=predictions,RE=pred_effects,FE=fix_effects)
  
  
  
  
  fix_train_PC <- as.matrix(cbind(myCV_train,myPCA_train))
  fix_test_PC  <- as.matrix(cbind(myCV_test,myPCA_test))
  p <- ncol(fix_train_PC)
  XtX <- crossprod(fix_train_PC, fix_train_PC)
  rank.X <- qr(XtX)$rank
  if (rank.X < p) {
    sm2=findLinearCombos(fix_train_PC)$remove
    fix_train_PC=fix_train_PC[,-sm2]
    fix_test_PC=fix_test_PC[,-sm2]}
  gc()
  rrBLUP_model_PC <- mixed.solve(y = myY_train,
                                 Z = myGD_train,
                                 X = fix_train_PC)
  
  pred_effects_PC <- myGD_test %*% rrBLUP_model_PC$u
  fix_effects_PC <- fix_test_PC  %*% rrBLUP_model_PC$beta
  predictions_PC <- c(pred_effects_PC) + c(fix_effects_PC)
  acc_PC <- cor(predictions_PC, myY_test, use = "pairwise.complete")
  sacc_PC <- cor(predictions_PC, myY_test, use = "pairwise.complete", method = c("spearman"))
  metrics_PC=postResample(pred=predictions_PC,obs=myY_test)
  results_PC=c(ACC_PC=acc_PC,SACC_PC=sacc_PC,metrics)
  prediction_PC=data.frame(test_phenotype,GEBV=predictions_PC,RE=pred_effects_PC,FE=fix_effects_PC)
  Accuracy=c(results,results_PC)
  names(Accuracy)<-c("Pearson","Spearman","RMSE","R2","MAE", "Pearson_PC","Spearman_PC","RMSE_PC","R2_PC","MAE_PC")
  Predictions<-cbind(prediction,prediction_PC[,3:5])
  names(Predictions)[6:8]<-c("GEBV_PC","RE_PC","FE_PC")
  results_ALL=list(Accuracy,Predictions)
  return(results_ALL)
}

MAS_vs_pc_mean_GAGS <- function(train_genotypes, train_phenotype,train_Y=NULL,train_GM=NULL,train_GD=NULL,train_PCA=NULL,test_genotypes,test_phenotype,test_PCA=NULL,GWAS="BLINK",alpha=0.05, folds = 5)
{
  library(rrBLUP)
  library(tidyr)
  library(Hmisc)
  library(caret)
  
  myY_train <- train_phenotype[,2]
  myY_test <- test_phenotype[,2]
  
  
  phenotype <- c(train_phenotype[,2],test_phenotype[,2])
  pheno_train <- c(train=train_phenotype[,2],test=test_phenotype[,2])
  pheno_train[grep("test",names(pheno_train),value=FALSE)] <- NA
  genotypes<-rbind(train_genotypes,test_genotypes)
  PCA<-rbind(train_PCA,test_PCA)
  
  GWASR<- GAPIT(Y = train_Y,
                GD = train_GD,
                GM = train_GM,
                PCA.total=3,
                model = GWAS,
                file.output=F)
  
  GWASSM <- which(GWASR$GWAS$P.value <= alpha/length(GWASR$GWAS$P.value))
  if(length(GWASSM)==0){
   
    
    acc <- NA
    sacc <- NA
    metrics=c(RMSE=NA,Rsquared=NA,MAE=NA)
    results=c(ACC=acc,SACC=sacc,metrics)
    prediction=NA
    
    acc_PC <- NA
    sacc_PC <- NA
    metrics_PC=c(RMSE=NA,Rsquared=NA,MAE=NA)
    results_PC=c(ACC_PC=acc_PC,SACC_PC=sacc_PC,metrics)
    
    prediction_PC=NA
    
    Accuracy=c(results,results_PC)
    
    names(Accuracy)<-c("Pearson","Spearman","RMSE","R2","MAE", "Pearson_PC","Spearman_PC","RMSE_PC","R2_PC","MAE_PC")
    Predictions<-cbind(prediction,prediction_PC)
    
    results_ALL=list(Accuracy,Predictions)
    
  }else{
    sm=as.character(GWASR$GWAS[GWASSM,]$SNP)
    
    CV <- genotypes[,sm]
    

    MAS_train <- data.frame(CV)
    MAS_test  <- data.frame(CV)
    
    GLM_data <- data.frame(pheno_train, MAS_train)

    names(GLM_data)[1] <- "Y"
    #Linear model to calculate effects
    #You can run all signficant markers at once to see cumulative effect
    MAS_model <- lm(Y ~ ., data = GLM_data)
    MAS_pred <- predict(MAS_model, MAS_test)
    
    acc <- cor(MAS_pred[-c(1:length(train_phenotype[,2]))], myY_test, use = "pairwise.complete")
    sacc <- cor(MAS_pred[-c(1:length(train_phenotype[,2]))], myY_test, use = "pairwise.complete", method = c("spearman"))
    metrics=postResample(pred=MAS_pred[-c(1:length(train_phenotype[,2]))],obs=myY_test)
    results=c(ACC=acc,SACC=sacc,metrics)
    prediction=data.frame(test_phenotype,GEBV=MAS_pred[-c(1:length(train_phenotype[,2]))])
    
    
    MAS_train_PC <- data.frame(cbind(CV,PCA))
    MAS_test_PC  <- data.frame(cbind(CV,PCA))
    
    GLM_data_PC <- data.frame(pheno_train, MAS_train_PC)

    names(GLM_data_PC)[1] <- "Y"
    #Linear model to calculate effects
    #You can run all signficant markers at once to see cumulative effect
    MAS_model_PC <- lm(Y ~ ., data = GLM_data_PC)
    MAS_pred_PC <- predict(MAS_model_PC, MAS_test_PC)
    
    acc_PC <- cor(MAS_pred_PC[-c(1:length(train_phenotype[,2]))], myY_test, use = "pairwise.complete")
    sacc_PC <- cor(MAS_pred_PC[-c(1:length(train_phenotype[,2]))], myY_test, use = "pairwise.complete", method = c("spearman"))
    metrics_PC=postResample(pred=MAS_pred_PC[-c(1:length(train_phenotype[,2]))],obs=myY_test)
    results_PC=c(ACC_PC=acc_PC,SACC_PC=sacc_PC,metrics)
    prediction_PC=data.frame(test_phenotype,GEBV=MAS_pred_PC[-c(1:length(train_phenotype[,2]))])
    
    Accuracy=c(results,results_PC)
    names(Accuracy)<-c("Pearson","Spearman","RMSE","R2","MAE", "Pearson_PC","Spearman_PC","RMSE_PC","R2_PC","MAE_PC")
    Predictions<-cbind(prediction,prediction_PC[,3])
    names(Predictions)[4]<-c("GEBV_PC")
    results_ALL=list(Accuracy,Predictions)
  }
  
  return(results_ALL)
}

MAS_vs_pc_mean_GAGS_10 <- function(train_genotypes, train_phenotype,train_Y=NULL,train_GM=NULL,train_GD=NULL,train_PCA=NULL,test_genotypes,test_phenotype,test_PCA=NULL,GWAS="BLINK",alpha=0.05, folds = 5)
{
  library(rrBLUP)
  library(tidyr)
  library(Hmisc)
  library(caret)
  
  
  myY_train <- train_phenotype[,2]
  myY_test <- test_phenotype[,2]
  
  
  phenotype <- c(train_phenotype[,2],test_phenotype[,2])
  pheno_train <- c(train=train_phenotype[,2],test=test_phenotype[,2])
  pheno_train[grep("test",names(pheno_train),value=FALSE)] <- NA
  genotypes<-rbind(train_genotypes,test_genotypes)
  PCA<-rbind(train_PCA,test_PCA)
  
  GWASR<- GAPIT(Y = train_Y,
                GD = train_GD,
                GM = train_GM,
                PCA.total=3,
                model = GWAS,
                file.output=F)
  
  top10=GWASR$GWAS[order(GWASR$GWAS$P.value),]
  GWASSM=top10[1:10,]$SNP

  
  CV <- genotypes[,GWASSM]
  
  
  MAS_train <- data.frame(CV)
  MAS_test  <- data.frame(CV)
  
  
  GLM_data <- data.frame(pheno_train, MAS_train)

  names(GLM_data)[1] <- "Y"
  #Linear model to calculate effects
  #You can run all signficant markers at once to see cumulative effect
  MAS_model <- lm(Y ~ ., data = GLM_data)
  MAS_pred <- predict(MAS_model, MAS_test)
  
  acc <- cor(MAS_pred[-c(1:length(train_phenotype[,2]))], myY_test, use = "pairwise.complete")
  sacc <- cor(MAS_pred[-c(1:length(train_phenotype[,2]))], myY_test, use = "pairwise.complete", method = c("spearman"))
  metrics=postResample(pred=MAS_pred[-c(1:length(train_phenotype[,2]))],obs=myY_test)
  results=c(ACC=acc,SACC=sacc,metrics)
  prediction=data.frame(test_phenotype,GEBV=MAS_pred[-c(1:length(train_phenotype[,2]))])
  
  
  MAS_train_PC <- data.frame(cbind(CV,PCA))
  MAS_test_PC  <- data.frame(cbind(CV,PCA))
  
  GLM_data_PC <- data.frame(pheno_train, MAS_train_PC)

  names(GLM_data_PC)[1] <- "Y"
  #Linear model to calculate effects
  #You can run all signficant markers at once to see cumulative effect
  MAS_model_PC <- lm(Y ~ ., data = GLM_data_PC)
  MAS_pred_PC <- predict(MAS_model_PC, MAS_test_PC)
  
  acc_PC <- cor(MAS_pred_PC[-c(1:length(train_phenotype[,2]))], myY_test, use = "pairwise.complete")
  sacc_PC <- cor(MAS_pred_PC[-c(1:length(train_phenotype[,2]))], myY_test, use = "pairwise.complete", method = c("spearman"))
  metrics_PC=postResample(pred=MAS_pred_PC[-c(1:length(train_phenotype[,2]))],obs=myY_test)
  results_PC=c(ACC_PC=acc_PC,SACC_PC=sacc_PC,metrics)
  prediction_PC=data.frame(test_phenotype,GEBV=MAS_pred_PC[-c(1:length(train_phenotype[,2]))])
  
  Accuracy=c(results,results_PC)
  names(Accuracy)<-c("Pearson","Spearman","RMSE","R2","MAE", "Pearson_PC","Spearman_PC","RMSE_PC","R2_PC","MAE_PC")
  Predictions<-cbind(prediction,prediction_PC[,3])
  names(Predictions)[4]<-c("GEBV_PC")
  results_ALL=list(Accuracy,Predictions)
  return(results_ALL)
}


MAS_vs_pc_mean <- function(train_phenotype,train_PCA=NULL,train_CV=NULL,test_phenotype,test_PCA=NULL,test_CV=NULL,GWAS="BLINK",alpha=0.05, folds = 5)
{
  library(rrBLUP)
  library(tidyr)
  library(Hmisc)
  library(caret)
  
  myY_train <- train_phenotype[,2]
  myY_test <- test_phenotype[,2]
  phenotype <- c(train_phenotype[,2],test_phenotype[,2])
  pheno_train <- c(train=train_phenotype[,2],test=test_phenotype[,2])
  pheno_train[grep("test",names(pheno_train),value=FALSE)] <- NA
  if(length(ncol(train_CV))==0){
    CV<-c(train_CV,test_CV)
  }else{
    CV<-rbind(train_CV,test_CV)
  }
  
  PCA<-rbind(train_PCA,test_PCA)
  
  
  MAS_train <- data.frame(CV)
  MAS_test  <- data.frame(CV)
  
  GLM_data <- data.frame(pheno_train, MAS_train)

  names(GLM_data)[1] <- "Y"
  #Linear model to calculate effects
  #You can run all signficant markers at once to see cumulative effect
  MAS_model <- lm(Y ~ ., data = GLM_data)
  MAS_pred <- predict(MAS_model, MAS_test)
  
  acc <- cor(MAS_pred[-c(1:length(train_phenotype[,2]))], myY_test, use = "pairwise.complete")
  sacc <- cor(MAS_pred[-c(1:length(train_phenotype[,2]))], myY_test, use = "pairwise.complete", method = c("spearman"))
  metrics=postResample(pred=MAS_pred[-c(1:length(train_phenotype[,2]))],obs=myY_test)
  results=c(ACC=acc,SACC=sacc,metrics)
  prediction=data.frame(test_phenotype,GEBV=MAS_pred[-c(1:length(train_phenotype[,2]))])
  
  
  MAS_train_PC <- data.frame(cbind(CV,PCA))
  MAS_test_PC  <- data.frame(cbind(CV,PCA))
  
  GLM_data_PC <- data.frame(pheno_train, MAS_train_PC)

  names(GLM_data_PC)[1] <- "Y"
  #Linear model to calculate effects
  #You can run all signficant markers at once to see cumulative effect
  MAS_model_PC <- lm(Y ~ ., data = GLM_data_PC)
  MAS_pred_PC <- predict(MAS_model_PC, MAS_test_PC)
  
  acc_PC <- cor(MAS_pred_PC[-c(1:length(train_phenotype[,2]))], myY_test, use = "pairwise.complete")
  sacc_PC <- cor(MAS_pred_PC[-c(1:length(train_phenotype[,2]))], myY_test, use = "pairwise.complete", method = c("spearman"))
  metrics_PC=postResample(pred=MAS_pred_PC[-c(1:length(train_phenotype[,2]))],obs=myY_test)
  results_PC=c(ACC_PC=acc_PC,SACC_PC=sacc_PC,metrics)
  prediction_PC=data.frame(test_phenotype,GEBV=MAS_pred_PC[-c(1:length(train_phenotype[,2]))])
  
  Accuracy=c(results,results_PC)
  names(Accuracy)<-c("Pearson","Spearman","RMSE","R2","MAE", "Pearson_PC","Spearman_PC","RMSE_PC","R2_PC","MAE_PC")
  Predictions<-cbind(prediction,prediction_PC[,3])
  names(Predictions)[4]<-c("GEBV_PC")
  results_ALL=list(Accuracy,Predictions)
  return(results_ALL)
}

BLUP_pc_mean_GAGS_QTN <- function(genotypes, phenotype,Y=NULL,GM=NULL,GD=NULL,PCA=NULL,CV=NULL,GWAS="BLINK",alpha=0.05, folds = 5,QTN=10)
{
  library(rrBLUP)
  library(tidyr)
  library(Hmisc)
  library(caret)
  # Make the CV list
  fold_list <- make_CV_sets(length(phenotype[,2]), k = folds)
  
  BGLR_acc_results <- list()
  BGLR_acc_predictions<- list()
  Predictions_ALL<-list()
  for (i in 1:length(fold_list)){
    fold_indices <- which(fold_list[[i]])
    
    # Split into training and testing data
    # Calculate the GS model using rrBLUP
    myGD_train <- as.matrix(genotypes[fold_indices,])
    myGD_train=apply(myGD_train,2,function(x) recode(x,"0"="-1","1"="0","2"="1"))
    myGD_train=apply(myGD_train,2,as.numeric)
    myGD_test <- as.matrix(genotypes[-fold_indices,])
    myGD_test=apply(myGD_test,2,function(x) recode(x,"0"="-1","1"="0","2"="1"))
    myGD_test=apply(myGD_test,2,as.numeric)
    
    GD_train <- GD[fold_indices,]
    GD_test <- GD[-fold_indices,]
    Y_train <-Y[fold_indices,c(1,2)]
    Y_test <-Y[-fold_indices,c(1,2)]
    myY_train <- phenotype[fold_indices,2]
    myY_test <- phenotype[-fold_indices,2]
    
    GWASR<- GAPIT(Y = Y_train,
                  GD = GD_train,
                  GM = GM,
                  PCA.total=3,
                  model = GWAS,
                  file.output=F)
    
    top10=GWASR$GWAS[order(GWASR$GWAS$P.value),]
    GWASSM=top10[1:QTN,]$SNP
    myCV_train <- myGD_train[,GWASSM]
    myCV_test <- myGD_test[,GWASSM]
    myPCA_train <- PCA[fold_indices,]
    myPCA_test <- PCA[-fold_indices,]
    
    fix_train <- as.matrix(myCV_train)
    fix_test  <- as.matrix(myCV_test)
    p <- ncol(fix_train)
    XtX <- crossprod(fix_train, fix_train)
    rank.X <- qr(XtX)$rank
    if (rank.X < p) {
      sm2=findLinearCombos(fix_train)$remove
      fix_train=fix_train[,-sm2]
      fix_test=fix_test[,-sm2]}
    gc()
    rrBLUP_model <- mixed.solve(y = myY_train,
                                Z = myGD_train,
                                X = fix_train)
    
    pred_effects <- myGD_test %*% rrBLUP_model$u
    fix_effects <- fix_test  %*% rrBLUP_model$beta
    predictions <- c(pred_effects) + c(fix_effects)
    acc <- cor(predictions, myY_test, use = "pairwise.complete")
    sacc <- cor(predictions, myY_test, use = "pairwise.complete", method = c("spearman"))
    metrics=postResample(pred=predictions,obs=myY_test)
    results=c(ACC=acc,SACC=sacc,metrics)
    prediction=data.frame(Y_test,GEBV=predictions,RE=pred_effects,FE=fix_effects)
    
    
    
    
    fix_train_PC <- as.matrix(cbind(myCV_train,myPCA_train))
    fix_test_PC  <- as.matrix(cbind(myCV_test,myPCA_test))
    p <- ncol(fix_train_PC)
    XtX <- crossprod(fix_train_PC, fix_train_PC)
    rank.X <- qr(XtX)$rank
    if (rank.X < p) {
      sm2=findLinearCombos(fix_train_PC)$remove
      fix_train_PC=fix_train_PC[,-sm2]
      fix_test_PC=fix_test_PC[,-sm2]}
    gc()
    rrBLUP_model_PC <- mixed.solve(y = myY_train,
                                   Z = myGD_train,
                                   X = fix_train_PC)
    
    pred_effects_PC <- myGD_test %*% rrBLUP_model_PC$u
    fix_effects_PC <- fix_test_PC  %*% rrBLUP_model_PC$beta
    predictions_PC <- c(pred_effects_PC) + c(fix_effects_PC)
    acc_PC <- cor(predictions_PC, myY_test, use = "pairwise.complete")
    sacc_PC <- cor(predictions_PC, myY_test, use = "pairwise.complete", method = c("spearman"))
    metrics_PC=postResample(pred=predictions_PC,obs=myY_test)
    results_PC=c(ACC_PC=acc_PC,SACC_PC=sacc_PC,metrics)
    prediction_PC=data.frame(Y_test,GEBV=predictions_PC,RE=pred_effects_PC,FE=fix_effects_PC)
    
    
    
    
    Predictions<-cbind(prediction,prediction_PC[,3:5])
    Predictions_ALL[[i]]=Predictions
    BGLR_acc_results[[i]] <- list(results,results_PC)
    #BGLR_acc_predictions[[i]] <- list(results,results_PC,prediction,prediction_PC)
  }
  #, GBLUP_acc, "GBLUP"
  model_vect <- c("Pearson","Spearman","RMSE","R2","MAE", "Pearson_PC","Spearman_PC","RMSE_PC","R2_PC","MAE_PC")
  BGLR_acc_table <- data.frame(matrix(nrow = 0, ncol = 3))
  for (i in 1:length(BGLR_acc_results))
  {
    results_long <- data.frame(rep(i, length(model_vect)), model_vect, unlist(BGLR_acc_results[[i]]))
    BGLR_acc_table <- rbind(BGLR_acc_table, results_long)
  }
  names(BGLR_acc_table) <- c("fold", "model", "r")
  data_wide <- spread(BGLR_acc_table, model, r)
  acc_fold=data.frame(data_wide)
  results=colMeans(acc_fold[2:11], na.rm=TRUE)
  names(Predictions)[6:8]<-c("GEBV_PC","RE_PC","FE_PC")
  results_ALL=list(results,Predictions_ALL)
  return(results_ALL)
}

MAS_GWAS_QTN_CV_Mean <- function(genotypes, phenotype,Y=NULL,GM=NULL,GD=NULL,PCA=NULL,CV=NULL,GWAS="BLINK",alpha=0.05, folds = 5,QTN=10)
{
  library(rrBLUP)
  library(tidyr)
  library(Hmisc)
  library(caret)
  library(dplyr)
  # Make the CV list
  fold_list <- make_CV_sets(length(phenotype[,2]), k = folds)
  
  BGLR_acc_results <- list()
  BGLR_acc_predictions<- list()
  Predictions_ALL<-list()
  for (i in 1:length(fold_list)){
    fold_indices <- which(fold_list[[i]])
    
    # Split into training and testing data
    # Calculate the GS model using rrBLUP
    myGD_train <- as.matrix(genotypes[fold_indices,])
    myGD_train=apply(myGD_train,2,function(x) recode(x,"0"="-1","1"="0","2"="1"))
    myGD_train=apply(myGD_train,2,as.numeric)
    myGD_test <- as.matrix(genotypes[-fold_indices,])
    myGD_test=apply(myGD_test,2,function(x) recode(x,"0"="-1","1"="0","2"="1"))
    myGD_test=apply(myGD_test,2,as.numeric)
    
    GD_train <- GD[fold_indices,]
    GD_test <- GD[-fold_indices,]
    Y_train <-Y[fold_indices,c(1,2)]
    Y_test <-Y[-fold_indices,c(1,2)]
    
    pheno_train <- phenotype[,2]
    pheno_train[-fold_indices] <- NA
    
    myY_train <- phenotype[fold_indices,2]
    myY_test <- phenotype[-fold_indices,2]
    
    GWASR<- GAPIT(Y = Y_train,
                  GD = GD_train,
                  GM = GM,
                  PCA.total=3,
                  model = GWAS,
                  file.output=F)
    
    top10=GWASR$GWAS[order(GWASR$GWAS$P.value),]
    GWASSM=top10[1:QTN,]$SNP
    
    
    CV <- genotypes[,GWASSM]
    myPCA_train <- PCA
    myPCA_test <- PCA
    
    
    MAS_train <- data.frame(CV)
    MAS_test  <- data.frame(CV)
    
    GLM_data <- data.frame(pheno_train, MAS_train)
    
    names(GLM_data)[1] <- "Y"
    #Linear model to calculate effects
    #You can run all signficant markers at once to see cumulative effect
    MAS_model <- lm(Y ~ ., data = GLM_data)
    MAS_pred <- predict(MAS_model, MAS_test)
    
    acc <- cor(MAS_pred[-fold_indices], myY_test, use = "pairwise.complete")
    sacc <- cor(MAS_pred[-fold_indices], myY_test, use = "pairwise.complete", method = c("spearman"))
    metrics=postResample(pred=MAS_pred[-fold_indices],obs=myY_test)
    results=c(ACC=acc,SACC=sacc,metrics)
    prediction=data.frame(Y_test,GEBV=MAS_pred[-fold_indices])
    
    
    
    
    MAS_train_PC <- data.frame(cbind(MAS_train,myPCA_train))
    MAS_test_PC  <- data.frame(cbind(MAS_test,myPCA_test))
    
    GLM_data_PC <- data.frame(pheno_train, MAS_train_PC)
    
    names(GLM_data_PC)[1] <- "Y"
    #Linear model to calculate effects
    #You can run all signficant markers at once to see cumulative effect
    MAS_model_PC <- lm(Y ~ ., data = GLM_data_PC)
    MAS_pred_PC <- predict(MAS_model_PC, MAS_test_PC)
    
    acc_PC <- cor(MAS_pred_PC[-fold_indices], myY_test, use = "pairwise.complete")
    sacc_PC <- cor(MAS_pred_PC[-fold_indices], myY_test, use = "pairwise.complete", method = c("spearman"))
    metrics_PC=postResample(pred=MAS_pred_PC[-fold_indices],obs=myY_test)
    results_PC=c(ACC_PC=acc_PC,SACC_PC=sacc_PC,metrics)
    prediction_PC=data.frame(Y_test,GEBV=MAS_pred_PC[-fold_indices])
    
    
    Predictions<-cbind(prediction,prediction_PC[,3])
    Predictions_ALL[[i]]=Predictions
    BGLR_acc_results[[i]] <- list(results,results_PC)
    #BGLR_acc_predictions[[i]] <- list(results,results_PC,prediction,prediction_PC)
  }
  #, GBLUP_acc, "GBLUP"
  model_vect <- c("Pearson","Spearman","RMSE","R2","MAE", "Pearson_PC","Spearman_PC","RMSE_PC","R2_PC","MAE_PC")
  BGLR_acc_table <- data.frame(matrix(nrow = 0, ncol = 3))
  for (i in 1:length(BGLR_acc_results))
  {
    results_long <- data.frame(rep(i, length(model_vect)), model_vect, unlist(BGLR_acc_results[[i]]))
    BGLR_acc_table <- rbind(BGLR_acc_table, results_long)
  }
  names(BGLR_acc_table) <- c("fold", "model", "r")
  data_wide <- spread(BGLR_acc_table, model, r)
  acc_fold=data.frame(data_wide)
  results=colMeans(acc_fold[2:11], na.rm=TRUE)
  names(Predictions)[4]<-c("GEBV_PC")
  results_ALL=list(results,Predictions_ALL)
  return(results_ALL)
}

BLUP_vs_pc_mean_GAGS_QTN <- function(train_genotypes, train_phenotype,train_Y=NULL,train_GM=NULL,train_GD=NULL,train_PCA=NULL,test_genotypes,test_phenotype,test_PCA=NULL,GWAS="BLINK",alpha=0.05, folds = 5,QTN=10)
{
  library(rrBLUP)
  library(tidyr)
  library(Hmisc)
  library(caret)
  
  myGD_train <- as.matrix(train_genotypes)
  myGD_train=apply(myGD_train,2,function(x) recode(x,"0"="-1","1"="0","2"="1"))
  myGD_train=apply(myGD_train,2,as.numeric)
  myGD_test <- as.matrix(test_genotypes)
  myGD_test=apply(myGD_test,2,function(x) recode(x,"0"="-1","1"="0","2"="1"))
  myGD_test=apply(myGD_test,2,as.numeric)
  
  myY_train <- train_phenotype[,2]
  myY_test <- test_phenotype[,2]
  
  GWASR<- GAPIT(Y = train_Y,
                GD = train_GD,
                GM = train_GM,
                PCA.total=3,
                model = GWAS,
                file.output=F)
  
  top10=GWASR$GWAS[order(GWASR$GWAS$P.value),]
  GWASSM=top10[1:QTN,]$SNP
  
  
  myCV_train <- myGD_train[,GWASSM]
  myCV_test <- myGD_test[,GWASSM]
  myPCA_train <- train_PCA
  myPCA_test <- test_PCA
  
  
  fix_train <- as.matrix(myCV_train)
  fix_test  <- as.matrix(myCV_test)
  p <- ncol(fix_train)
  XtX <- crossprod(fix_train, fix_train)
  rank.X <- qr(XtX)$rank
  if (rank.X < p) {
    sm2=findLinearCombos(fix_train)$remove
    fix_train=fix_train[,-sm2]
    fix_test=fix_test[,-sm2]}
  gc()
  rrBLUP_model <- mixed.solve(y = myY_train,
                              Z = myGD_train,
                              X = fix_train)
  
  pred_effects <- myGD_test %*% rrBLUP_model$u
  fix_effects <- fix_test  %*% rrBLUP_model$beta
  predictions <- c(pred_effects) + c(fix_effects)
  acc <- cor(predictions, myY_test, use = "pairwise.complete")
  sacc <- cor(predictions, myY_test, use = "pairwise.complete", method = c("spearman"))
  metrics=postResample(pred=predictions,obs=myY_test)
  results=c(ACC=acc,SACC=sacc,metrics)
  prediction=data.frame(test_phenotype,GEBV=predictions,RE=pred_effects,FE=fix_effects)
  
  
  
  
  fix_train_PC <- as.matrix(cbind(myCV_train,myPCA_train))
  fix_test_PC  <- as.matrix(cbind(myCV_test,myPCA_test))
  p <- ncol(fix_train_PC)
  XtX <- crossprod(fix_train_PC, fix_train_PC)
  rank.X <- qr(XtX)$rank
  if (rank.X < p) {
    sm2=findLinearCombos(fix_train_PC)$remove
    fix_train_PC=fix_train_PC[,-sm2]
    fix_test_PC=fix_test_PC[,-sm2]}
  gc()
  rrBLUP_model_PC <- mixed.solve(y = myY_train,
                                 Z = myGD_train,
                                 X = fix_train_PC)
  
  pred_effects_PC <- myGD_test %*% rrBLUP_model_PC$u
  fix_effects_PC <- fix_test_PC  %*% rrBLUP_model_PC$beta
  predictions_PC <- c(pred_effects_PC) + c(fix_effects_PC)
  acc_PC <- cor(predictions_PC, myY_test, use = "pairwise.complete")
  sacc_PC <- cor(predictions_PC, myY_test, use = "pairwise.complete", method = c("spearman"))
  metrics_PC=postResample(pred=predictions_PC,obs=myY_test)
  results_PC=c(ACC_PC=acc_PC,SACC_PC=sacc_PC,metrics)
  prediction_PC=data.frame(test_phenotype,GEBV=predictions_PC,RE=pred_effects_PC,FE=fix_effects_PC)
  Accuracy=c(results,results_PC)
  names(Accuracy)<-c("Pearson","Spearman","RMSE","R2","MAE", "Pearson_PC","Spearman_PC","RMSE_PC","R2_PC","MAE_PC")
  Predictions<-cbind(prediction,prediction_PC[,3:5])
  names(Predictions)[6:8]<-c("GEBV_PC","RE_PC","FE_PC")
  results_ALL=list(Accuracy,Predictions)
  return(results_ALL)
}

MAS_vs_pc_mean_GAGS_QTN <- function(train_genotypes, train_phenotype,train_Y=NULL,train_GM=NULL,train_GD=NULL,train_PCA=NULL,test_genotypes,test_phenotype,test_PCA=NULL,GWAS="BLINK",alpha=0.05, folds = 5,QTN=10)
{
  library(rrBLUP)
  library(tidyr)
  library(Hmisc)
  library(caret)
  
  
  myY_train <- train_phenotype[,2]
  myY_test <- test_phenotype[,2]
  
  
  phenotype <- c(train_phenotype[,2],test_phenotype[,2])
  pheno_train <- c(train=train_phenotype[,2],test=test_phenotype[,2])
  pheno_train[grep("test",names(pheno_train),value=FALSE)] <- NA
  genotypes<-rbind(train_genotypes,test_genotypes)
  PCA<-rbind(train_PCA,test_PCA)
  
  GWASR<- GAPIT(Y = train_Y,
                GD = train_GD,
                GM = train_GM,
                PCA.total=3,
                model = GWAS,
                file.output=F)
  
  top10=GWASR$GWAS[order(GWASR$GWAS$P.value),]
  GWASSM=top10[1:QTN,]$SNP
  
  
  CV <- genotypes[,GWASSM]
  
  
  MAS_train <- data.frame(CV)
  MAS_test  <- data.frame(CV)
  
  
  GLM_data <- data.frame(pheno_train, MAS_train)
  
  names(GLM_data)[1] <- "Y"
  #Linear model to calculate effects
  #You can run all signficant markers at once to see cumulative effect
  MAS_model <- lm(Y ~ ., data = GLM_data)
  MAS_pred <- predict(MAS_model, MAS_test)
  
  acc <- cor(MAS_pred[-c(1:length(train_phenotype[,2]))], myY_test, use = "pairwise.complete")
  sacc <- cor(MAS_pred[-c(1:length(train_phenotype[,2]))], myY_test, use = "pairwise.complete", method = c("spearman"))
  metrics=postResample(pred=MAS_pred[-c(1:length(train_phenotype[,2]))],obs=myY_test)
  results=c(ACC=acc,SACC=sacc,metrics)
  prediction=data.frame(test_phenotype,GEBV=MAS_pred[-c(1:length(train_phenotype[,2]))])
  
  
  MAS_train_PC <- data.frame(cbind(CV,PCA))
  MAS_test_PC  <- data.frame(cbind(CV,PCA))
  
  GLM_data_PC <- data.frame(pheno_train, MAS_train_PC)
  
  names(GLM_data_PC)[1] <- "Y"
  #Linear model to calculate effects
  #You can run all signficant markers at once to see cumulative effect
  MAS_model_PC <- lm(Y ~ ., data = GLM_data_PC)
  MAS_pred_PC <- predict(MAS_model_PC, MAS_test_PC)
  
  acc_PC <- cor(MAS_pred_PC[-c(1:length(train_phenotype[,2]))], myY_test, use = "pairwise.complete")
  sacc_PC <- cor(MAS_pred_PC[-c(1:length(train_phenotype[,2]))], myY_test, use = "pairwise.complete", method = c("spearman"))
  metrics_PC=postResample(pred=MAS_pred_PC[-c(1:length(train_phenotype[,2]))],obs=myY_test)
  results_PC=c(ACC_PC=acc_PC,SACC_PC=sacc_PC,metrics)
  prediction_PC=data.frame(test_phenotype,GEBV=MAS_pred_PC[-c(1:length(train_phenotype[,2]))])
  
  Accuracy=c(results,results_PC)
  names(Accuracy)<-c("Pearson","Spearman","RMSE","R2","MAE", "Pearson_PC","Spearman_PC","RMSE_PC","R2_PC","MAE_PC")
  Predictions<-cbind(prediction,prediction_PC[,3])
  names(Predictions)[4]<-c("GEBV_PC")
  results_ALL=list(Accuracy,Predictions)
  return(results_ALL)
}