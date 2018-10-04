#trials <- 10
maxsteps = 100
my.minlam <- 10^(-5)
nfold <- 10
outerCV_folds <- 5
print(paste('nfold=',nfold))
info <- vector(mode="character",length=3)
args <- commandArgs(trailingOnly = T)
#args <- c("DEEP_allDS_divergent_plus.txt","T")
#options(warn = -1)
a <- length(args)
params <- readLines(args[1])
if(F){
intercept_mode <- args[2]
if(intercept_mode == "FALSE" || intercept_mode == "F" || intercept_mode == "0"){
  intercept_mode <- F
}else if(intercept_mode == "TRUE" || intercept_mode == "T" || intercept_mode == "1"){
  intercept_mode <- T
}else{
  stop(paste("A TRUE or FALSE should be provided for the intercept_mode argument. The given argument by user is ",intercept_mode,sep=""))
}
ylab <- args[3]
}
ylab <- "measured"


print(params)
dataset.inf <- which(params=="Dataset (Provide the path to the datasets)")
gamma.inf <- which(params=="Model Optimization")
feature.inf <- which(params=="Feature Group Names")
output.inf <- which(params=="Output (Provide the path to the location where the results should be saved)")
if(dataset.inf > gamma.inf || gamma.inf > feature.inf || feature.inf > output.inf || dataset.inf > output.inf)
  stop("The argument file is in incorrect format. The order for the parameters must be like Datasets, Model Optimization, Feature Groups Names, and Output.")
dataset.tokens <- strsplit(params[seq(dataset.inf+1,gamma.inf-1,1)],split = '=')
inPaths <- vector('character',length=length(dataset.tokens))
for(i in seq(1,length(dataset.tokens),2)){
  if(!(gsub(" *$","",dataset.tokens[[i]][1])=="X" && gsub(" *$","",dataset.tokens[[(i+1)]][1])=="Y")){
    if(i==1)
      number.str <- "st"
    else if(i==2)
      number.str <- "nd"
    else
      number.str <- "th"
    stop(paste("Either X o Y or both is missing for the ",i,number.str," dataset!",sep=''))
  }
  inPaths[i] <- gsub(" *","",dataset.tokens[[i]][2])
  inPaths[(i+1)] <- gsub(" *","",dataset.tokens[[(i+1)]][2])
}

gamma <- list(start=NULL,end=NULL,step=NULL)

gamma.tokens <- strsplit(params[seq(gamma.inf+1,feature.inf-1,1)],split = '=')
print(gamma.tokens[[1]][1])
print(gamma.tokens[[2]][1])
print(gamma.tokens[[3]][1])
if(length(gamma.tokens)==3){
  if(tolower(gsub(" * ","",gamma.tokens[[1]][1]))!='percent'){
    if(tolower(gsub(" * ","",gamma.tokens[[1]][1]))!='nfold'){
      if(tolower(gsub(" * ","",gamma.tokens[[1]][1]))=='interceptmode'){
        percent <- .8
        intercept_mode <- gsub(" * ","",gamma.tokens[[1]][2])

        if(intercept_mode == "FALSE" || intercept_mode == "F" || intercept_mode == "0"){
          intercept_mode <- F
        }else if(intercept_mode == "TRUE" || intercept_mode == "T" || intercept_mode == "1"){
          intercept_mode <- T
        }else{
          stop(paste("A TRUE or FALSE should be provided for the intercept_mode argument. The given argument by user is ",intercept_mode,sep=""))
        }
        nfold <- 10
        print('Warning! No percentage and nfold values are provided for model optimization. Setting percent to 0.8 and nfold to 10.')
      }else{
        percent <- 0.8
        nfold <- 10
        intercept_mode <- T
        print('Warning! No training optimization values are provided. Setting percent to 0.8, nfold to 10 and intercept mode to T')
      }
    }else{#nfold is specified as the first parameter for optimization
      if(tolower(gsub(" * ","",gamma.tokens[[2]][1])=='interceptmode')){
        percent <- 0.8
        nfold <- as.numeric(gamma.tokens[[1]][2])
        intercept_mode <- gsub(" * ","",gamma.tokens[[2]][2])
        if(intercept_mode == "FALSE" || intercept_mode == "F" || intercept_mode == "0"){
          intercept_mode <- F
        }else if(intercept_mode == "TRUE" || intercept_mode == "T" || intercept_mode == "1"){
          intercept_mode <- T
        }else{
          stop(paste("A TRUE or FALSE should be provided for the intercept_mode argument. The given argument by user is ",intercept_mode,sep=""))
        }
        print(paste('Warning! No value for percent is provided. Setting percent to 0.8',sep=''))
      }else{#step is not specified
        percent <- 0.8
        nfold <- as.numeric(gamma.tokens[[1]][2])
        intercept_mode <- T
        print('Warning! No percent and intercept mode values are provided. Setting percent to 0.8 and intercept mode to T')
      }
    }
  }else if(tolower(gsub(" * ","",gamma.tokens[[2]][1]))=='nfold' & tolower(gsub(" * ","",gamma.tokens[[3]][1]))=='interceptmode'){
    print("last if")
    percent <- as.numeric(gamma.tokens[[1]][2])
    nfold <- as.numeric(gamma.tokens[[2]][2])
    intercept_mode <- gsub(" * ","",gamma.tokens[[3]][2])
    if(intercept_mode == "FALSE" || intercept_mode == "F" || intercept_mode == "0"){
      intercept_mode <- F
    }else if(intercept_mode == "TRUE" || intercept_mode == "T" || intercept_mode == "1"){
      intercept_mode <- T
    }else{
      stop(paste("A TRUE or FALSE should be provided for the intercept_mode argument. The given argument by user is ",intercept_mode,sep=""))
    }

  }
}
feature.names <- params[seq(feature.inf+1,output.inf-1,1)]
outPath <- params[(output.inf+1)]
dir.create(outPath,showWarnings = F)
graphComponents <- length(feature.names)
print(graphComponents)
source('wrapper_functions_RSSmeasureModelSelection_gamma_parallel.R')
source('UtilityFunctions_sequentialCV.R')
datasets <- list()
fl.allModels <- list()
nl.allModels <- list()
#percent <- 0.8
bestgammas <-vector(mode='numeric',length=length(inPaths)/2)
bestGammasIdx <- 1
ctr <- 1
for(p in seq(1,length(inPaths),2)){
  set.seed(1)
  fileName <- unlist(strsplit(unlist(strsplit(inPaths[(p+1)],split = '/'))[length(unlist(strsplit(inPaths[(p+1)],split = '/')))],split='.y'))[1]
  if(p==1){
    datasets.names <- fileName
  }else{
    datasets.names <- c(datasets.names,fileName)
  }
  print(paste('Reading ',fileName,'...',sep=''))
  data.x <- log2(as.matrix(read.table(inPaths[p],header = F))+1)
  #data.y <- as.numeric(readLines(inPaths[(p+1)],))#[1:nrow(data.x)])
  data.y <- log2(as.numeric(readLines(inPaths[(p+1)],))+1)#[1:nrow(data.x)])
  shuffle <- sample(length(data.y))
  data.x <- data.x[shuffle,]
  data.y <- data.y[shuffle]
  print(class(data.y))
  print(paste('Done reading ',fileName,'...',sep=''))
  bin.cnt <- ncol(data.x)/graphComponents
  edgs <- NULL
  for(i in 1:graphComponents)
  {
    edgs.partial <- c(bin.cnt*(i-1)+1, rep((bin.cnt*(i-1)+2):(bin.cnt*i - 1),each=2),bin.cnt*i)
    #edg_i <- c(edgs.partial[seq(1,(length(edgs.partial)/2-2))],edgs.partial[seq(3+(length(edgs.partial)/2),length(edgs.partial))]) ## removes the NFR edges completely
    edg_i <- c(edgs.partial[seq(1,(length(edgs.partial)/2))],edgs.partial[seq(3+(length(edgs.partial)/2),length(edgs.partial))]) ## connects NFR to the upstream part only
    edgs <- c(edgs,edg_i)
    #edgs <- c(edgs,bin.cnt*(i-1)+1, rep((bin.cnt*(i-1)+2):(bin.cnt*i - 1),each=2),bin.cnt*i)
  }

  datasets[[ctr]] <- list(x=data.x,y=data.y)
  
  partition <- data.partition(data.x,data.y,percent=percent)
  if(!intercept_mode){
    #### Taken from fusedlasso paper:
    ## We also assume that the predictors are standardized to have mean 0
    ## and unit variance, and the outcome yi has mean 0. Hence we do
    ## not need an intercept in model (1).
    train.mean <- list(x=mean(partition$train$x),y=mean(partition$train$y))
    train.sd <- sd(partition$train$x)
    partition$train$x <- scale(partition$train$x)
    partition$train$y <- scale(partition$train$y,scale=F)
    partition$test$x <- (partition$test$x - train.mean$x)/train.sd
    partition$test$y <- (partition$test$y - train.mean$y)
  }else{
    #train.group.min <- sapply(seq(graphComponents),function(i)min(min(partition$train$x[,seq((i-1)*bin.cnt+1,i*bin.cnt)])))
    #train.group.max <- sapply(seq(graphComponents),function(i)max(max(partition$train$x[,seq((i-1)*bin.cnt+1,i*bin.cnt)])))
    train.group.rescaled <- group.rescale(partition$train$x,0,1,bin.cnt)
    partition$train$x <- train.group.rescaled$x
    partition$test$x <- group.rescale(partition$test$x,0,1,bin.cnt,train.group.rescaled$min,train.group.rescaled$max)$x
  }
  outer_CV_partitions <- get.outerCV.partitions(data.x,data.y,n.folds=outerCV_folds)
  print(gamma)
  print(dim(partition$train$x))
  gammas <- 0
  for(i in seq(-5,5)){
    gammas <- c(gammas,(10^i))
  }
#  gammas <- c(.5,1,2)

  #################################################################################
  ################ Run the outer CV to evaluate the unbiased error ################
  
  fl.outerCV <- sapply(seq(outerCV_folds),function(i)fusedlasso.main_ForLambdaInterpolation(outer_CV_partitions[[i]]$train$x,outer_CV_partitions[[i]]$train$y,bin.cnt,edgs,gammas,intercept_mode))
  fl <- fusedlasso.main_ForLambdaInterpolation(partition$train$x,partition$train$y,bin.cnt,edgs,gammas,intercept_mode)

  if(intercept_mode){
    edgs <- edgs + 1
    edgs <- c(1,1,edgs)
  }
  
  outer_CV_fl <- list()
  outer_CV_error <- NULL
  for(i in seq(outerCV_folds)){
    outerCV.best.gamma <- fl.outerCV[,i]$gamma.best
    if(intercept_mode){
      outer_CV_fl[[i]] <- fusedlasso(X=cbind(1,outer_CV_partitions[[i]]$train$x),y=outer_CV_partitions[[i]]$train$y,graph=graph(edgs,directed=F),gamma=outerCV.best.gamma,maxsteps = maxsteps,minlam=my.minlam)
    }else{
      outer_CV_fl[[i]] <- fusedlasso(X=outer_CV_partitions[[i]]$train$x,y=outer_CV_partitions[[i]]$train$y,graph=graph(edgs,directed=F),gamma=outerCV.best.gamma,maxsteps = maxsteps,minlam=my.minlam)
    }
    outerCV.min.lambda <- min(outer_CV_fl[[i]]$lambda) 
    beta.best <- coef.genlasso(outer_CV_fl[[i]],max(outerCV.min.lambda,fl$best.lambda))$beta
    if(intercept_mode){
      pred <- cbind(1,outer_CV_partitions[[i]]$train$x) %*% as.numeric(beta.best)
    }else{
      pred <- outer_CV_partitions[[i]]$train$x %*% as.numeric(beta.best)
    }
    err <- 1/length(pred) * sum((pred - outer_CV_partitions[[i]]$train$y)^2)
    outer_CV_error <- c(outer_CV_error,err)
  }
  writeLines(text=as.character(outer_CV_error),paste(outPath,'/outerCV_errors_',fileName,'.txt',sep=''))
  #################################################################################
  #################################################################################
  best.cv.gamma <- fl$gamma.best
  
  if(intercept_mode){
    fl.best <- fusedlasso(X=cbind(1,partition$train$x),y=partition$train$y,graph=graph(edgs,directed=F),gamma=best.cv.gamma,maxsteps = maxsteps,minlam=my.minlam)
  }else{
    fl.best <- fusedlasso(X=partition$train$x,y=partition$train$y,graph=graph(edgs,directed=F),gamma=best.cv.gamma,maxsteps = maxsteps,minlam=my.minlam)
  }
  fl.best.min.lambda <- min(fl.best$lambda) 
  beta.best <- coef.genlasso(fl.best,max(fl.best.min.lambda,fl$best.lambda))$beta
  print('Done running fusedlasso.main')
  cv.beta.mat.fl <- fl$cv.fl$cv.beta.mat
  print(c('dim(cv.beta.mat.fl)',dim(cv.beta.mat.fl)))
  cv.beta.mat.nl <- cv.beta.matrix.nl(partition$train$x,partition$train$y,nfold)
  print(c('save the fl, partition, and feature.names to RData files.'))
  save(fl,file=paste(outPath,'/fl_',fileName,'.RData',sep='')) 
  save(partition,file=paste(outPath,'/train_test_',fileName,'.RData',sep=''))
  save(feature.names,file=paste(outPath,'/featureNames_',fileName,'.RData',sep=''))

  fl.allModels[[ctr]] <- beta.best
  nl.allModels[[ctr]] <- cv.beta.mat.nl$best.nl

  bestgammas[bestGammasIdx] <- fl$gamma.best
  if(all(beta.best == 0)){
    print(bestgammas)
    stop('The trained Fussed lasso model contains all the coefficients as zero... The models is not trained properly (needs more parameter tuning, maybe!)')
  }

  save(cv.beta.mat.fl,cv.beta.mat.nl,file=paste(outPath,'/cv_beta_',fileName,'.RData',sep=''))
  pdf(paste(outPath,'stability_',fileName,'.pdf',sep=''))
  plot.stability(cv.beta.mat.fl,cv.beta.mat.nl$cv.beta,bin.cnt,feature.names,intercept_mode)
  plot.stability.var(cv.beta.mat.fl,cv.beta.mat.nl$cv.beta,intercept_mode)
  dev.off()

  ##########################################################################################
  print('preparing the details file')
  info[1] <- paste('best gamma(s)',paste(bestgammas,collapse='\t'))
  info[2] <- paste('median of fl',paste(fl.allModels[[ctr]],collapse="\t"))
  info[3] <- paste('median of nl',paste(nl.allModels[[ctr]],collapse="\t"))
  info[4] <- paste('number of CV folds',nfold,sep='\t')
  info[5] <- paste('number of maxsteps',maxsteps,sep='\t')
  writeLines(text=as.character(info),paste(outPath,'/details_',fileName,'.txt',sep=''))
  print('done writing the details file')
  ##########################################################################################
  
  correlations.fl <- matrix(0,ncol=3,nrow=2)#2 rows for spearman and pearson correlation and 3 columns for test, train, and entire dataset
  correlations.nl <- matrix(0,ncol=3,nrow=2)#2 rows for spearman and pearson correlation and 3 columns for test, train, and entire dataset
  
  ##########################################################################################
  #####################      plot coef heatmap      ########################################
  ##########################################################################################
  pdf(paste(outPath,'fusedLasso_coefs_heatmap',fileName,'.pdf',sep=''))
  print(fl$cv.fl$bestsol$beta)
  plot.histone.coef(beta.best,bin.cnt,feature.names,main=paste(fileName,'one model','gamma*',bestgammas[bestGammasIdx],sep='_'),intercept_mode)#,cluster_rows = F, cluster_cols = F)
  dev.off()

  pdf(paste(outPath,'standardLasso_coefs_heatmap',fileName,'.pdf',sep=''))
  plot.histone.coef(cv.beta.mat.nl$best.nl,bin.cnt,feature.names,main=paste(fileName,'one model','gamma*',bestgammas[bestGammasIdx],sep='_'),intercept_mode)#,cluster_rows = F, cluster_cols = F)
  dev.off()
  ##########################################################################################
  #####################      scatter plots test data      ##################################
  ##########################################################################################
  test.y <- partition$test$y
  if(intercept_mode){
    test.x <- partition$test$x
    train.x <- partition$train$x
    entire.x <- rbind(partition$train$x,partition$test$x)
  }else{
    test.x <- partition$test$x
    train.x <- partition$train$x
    entire.x <- rbind(train.x,test.x)
  }
  pdf(paste(outPath,'scatterPlots_test_',fileName,'.pdf',sep=''))
  plot.scatter(test.x,test.y,beta.best,0,ylab,main=paste('fl',fileName,sep='_'),is_fl=T,intercept_mode,xlab='prediction')
  ####################################################################

  plot.scatter(partition$test$x,partition$test$y,cv.beta.mat.nl$best.nl,0,ylab,main=paste('nl',fileName,sep='_'),is_fl=F,intercept_mode,xlab='prediction')
  dev.off()
  print('done test data scatter plots')
  ##########################################################################################
  #####################      scatter plots train data      ##################################
  ##########################################################################################
  pdf(paste(outPath,'scatterPlots_train_',fileName,'.pdf',sep=''))
  plot.scatter(train.x,partition$train$y,beta.best,0,ylab,main=paste('fl',fileName,sep='_'),is_fl=T,intercept_mode,xlab='prediction')
  ####################################################################

  plot.scatter(partition$train$x,partition$train$y,cv.beta.mat.nl$best.nl,0,ylab,main=paste('nl',fileName,sep='_'),is_fl=F,intercept_mode,xlab='prediction')
  dev.off()
  print('done train data scatter plots')
  ##########################################################################################
  #####################     scatter plots entire data      #################################
  ##########################################################################################
  pdf(paste(outPath,'scatterPlots_entireData_',fileName,'.pdf',sep=''))
  plot.scatter(entire.x,c(partition$train$y,partition$test$y),beta.best,0,ylab,main=paste('fl',fileName,sep='_'),is_fl=T,intercept_mode,xlab='prediction')

  ####################################################################

  plot.scatter(rbind(partition$train$x,partition$test$x),c(partition$train$y,partition$test$y),cv.beta.mat.nl$best.nl,0,ylab,main=paste('nl',fileName,sep='_'),is_fl=F,intercept_mode,xlab='prediction')
  dev.off()
  print('done entire data scatter plots')
  ##########################################################################################
  ##########################################################################################
  ##########################################################################################
  fl.beta <- fl$cv.fl$bestsol$beta
  if(intercept_mode)
    fl.beta <- fl$cv.fl$bestsol$beta[2:length(fl.beta)]
  beta.mat <- matrix(fl.beta,nrow=graphComponents,ncol=bin.cnt,byrow = T)
  pdf(paste(outPath,'barplots_',fileName,'.pdf',sep=''))
  print(coefficients.barplot(beta.mat,'','sum over the bins',feature.names))
  print(bins.barplot(beta.mat,'Bins','sum over features',seq(-floor(bin.cnt/2)+1,floor(bin.cnt/2),1)))
  dev.off()
  ctr <- ctr + 1
  bestGammasIdx <- bestGammasIdx + 1
}
save(datasets,file=paste(outPath,'/allDatasets.RData',sep=''))
save(fl.allModels,nl.allModels,file=paste(outPath,'/fl_nl_allModels.RData',sep=''))
if(length(datasets)>1){
  pdf(paste(outPath,'/cross_models.pdf',sep=''))
  plot.fl.accuracy.allmodels(fl.allModels,datasets,datasets.names,main='flasso test set',percent=percent,intercept_mode)
  plot.fl.accuracy.allmodels(nl.allModels,datasets,datasets.names,main='lasso test set',percent=percent,intercept_mode)
  dev.off()
}
print('best gamma(s)')
print(bestgammas)
print(warnings())
print(c('maxsteps',maxsteps))
