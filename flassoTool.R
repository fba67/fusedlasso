set.seed(0)
#trials <- 10
nfold <- 50
info <- vector(mode="character",length=3)
args <- commandArgs(trailingOnly = T)
#options(warn = -1)
a <- length(args)
params <- readLines(args[1])
ylab <- args[2]
print(params)
dataset.inf <- which(params=="Dataset (Provide the path to the datasets)")
gamma.inf <- which(params=="Gamma Range")
feature.inf <- which(params=="Feature Group Names")
output.inf <- which(params=="Output (Provide the path to the location where the results should be saved)")
if(dataset.inf > gamma.inf || gamma.inf > feature.inf || feature.inf > output.inf || dataset.inf > output.inf)
  stop("The argument file is in incorrect format. The order for the parameters must be like Datasets, Gamma Range, Feature Groups Names, and Output.")
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
#gamma <- list(start=0,end=1,step=.1)

gamma.tokens <- strsplit(params[seq(gamma.inf+1,feature.inf-1,1)],split = '=')
if(length(gamma.tokens)==3){
  if(tolower(gsub(" * ","",gamma.tokens[[1]][1]))!='start'){
    if(tolower(gsub(" * ","",gamma.tokens[[1]][1]))!='end'){
      if(tolower(gsub(" * ","",gamma.tokens[[1]][1]))=='step'){
        gamma$start <- 0
        gamma$step <- as.numeric(gamma.tokens[[1]][2])
        gamma$end <- 10*gamma$step
        print(paste('Warning! No Start and End values are provided for gamma. Setting Start to 0 and End to ',gamma$end,sep=''))
      }else{
        gamma$start <- 0
        gamma$step <- 0.1
        gamma$end <- 10*gamma$step
        print(paste('Warning! No Start and End, and Step values are provided for gamma. Setting Start to 0 and End to 1 and Step to 0.1',sep=''))
      }
    }else{#end is specified as the first parameter for gamma
      if(tolower(gsub(" * ","",gamma.tokens[[2]][1])=='step')){
        gamma$start <- 0
        gamma$step <- as.numeric(gamma.tokens[[2]][2])
        gamma$end <- as.numeric(gamma.tokens[[1]][2])
        print(paste('Warning! No Start value is provided for gamma. Setting Start to 0.',sep=''))
      }else{#step is not specified
        gamma$start <- 0
        gamma$end <- as.numeric(gamma.tokens[[1]][2])
        gamma$step <- gamma$end/10
        print(paste('Warning! No Start and Step values are provided for gamma. Setting Start to 0 and Step to ',gamma$step,sep=''))
      }
    }
  }else if(tolower(gsub(" * ","",gamma.tokens[[2]][1]))=='end' & tolower(gsub(" * ","",gamma.tokens[[3]][1]))=='step'){
    gamma$start <- as.numeric(gamma.tokens[[1]][2])
    gamma$end <- as.numeric(gamma.tokens[[2]][2])
    gamma$step <- as.numeric(gamma.tokens[[3]][2])
  }
}
feature.names <- params[seq(feature.inf+1,output.inf-1,1)]
outPath <- params[(output.inf+1)]
dir.create(outPath,showWarnings = F)
graphComponents <- length(feature.names)
print(graphComponents)
#source('wrapper_functions.R')
source('wrapper_functions_RSSmeasureModelSelection2.R')
datasets <- list()
fl.allModels <- list()
nl.allModels <- list()
percent <- 0.8
bestgammas <-vector(mode='numeric',length=length(inPaths)/2)
bestGammasIdx <- 1
ctr <- 1
for(p in seq(1,length(inPaths),2)){
  fileName <- unlist(strsplit(inPaths[p],split = '/'))[length(unlist(strsplit(inPaths[p],split = '/')))]
  if(p==1)
    datasets.names <- fileName
  else
    datasets.names <- c(datasets.names,fileName)
  print(paste('Reading ',fileName,'...',sep=''))
  data.x <- log2(as.matrix(read.table(inPaths[p],header = F))+1)
  #data.x <- scale(data.x)  
  #data.x <- data.x
  data.y <- as.numeric(readLines(inPaths[(p+1)],))#[1:nrow(data.x)])
  data.y <- log2(as.numeric(readLines(inPaths[(p+1)],))+1)#[1:nrow(data.x)])
  print(class(data.y))
  print(paste('Done reading ',fileName,'...',sep=''))
  bin.cnt <- ncol(data.x)/graphComponents
  edgs <- NULL
  for(i in 1:graphComponents)
  {
    edgs <- c(edgs,bin.cnt*(i-1)+1, rep((bin.cnt*(i-1)+2):(bin.cnt*i - 1),each=2),bin.cnt*i)
  }
  datasets[[ctr]] <- list(x=data.x,y=data.y)
  partition <- data.partition(data.x,data.y,percent=percent)
  print(gamma)
  print(dim(partition$train$x))
  fl <- fusedlasso.main(partition$train$x,partition$train$y,bin.cnt,edgs,seq(gamma$start,gamma$end,gamma$step))
  print('Done running fusedlasso.main')
  cv.beta.mat.fl <- fl$cv.fl$bestsol$cv.beta.mat
  print(c('dim(cv.beta.mat.fl)',dim(cv.beta.mat.fl)))
  cv.beta.mat.nl <- cv.beta.matrix.nl(partition$train$x,partition$train$y,nfold)
  save(fl,file=paste(outPath,'/fl_',fileName,'.RData',sep='')) 
  save(partition,file=paste(outPath,'/train_test_',fileName,'.RData',sep=''))
  save(feature.names,file=paste(outPath,'/featureNames_',fileName,'.RData',sep=''))

  fl.allModels[[ctr]] <- apply(cv.beta.mat.fl,2,FUN=median)
  nl.allModels[[ctr]] <- apply(cv.beta.mat.nl$cv.beta,2,FUN=median)

  bestgammas[bestGammasIdx] <- fl$gamma.best
  if(all(fl$cv.fl$bestsol$beta==0)){
    print(bestgammas)
    stop('The trained Fussed lasso model contains all the coefficients as zero... The models is not trained properly (needs more parameter tuning, maybe!)')
  }
  print(bestgammas)
  pdf(paste(outPath,'stability_',fileName,'.pdf',sep=''))

  plot.stability(cv.beta.mat.fl,cv.beta.mat.nl$cv.beta,bin.cnt,feature.names)
  plot.stability.var(cv.beta.mat.fl,cv.beta.mat.nl$cv.beta)
  dev.off()
  save(cv.beta.mat.fl,cv.beta.mat.nl,file=paste(outPath,'/cv_beta_',fileName,'.RData',sep=''))

  print('preparing the details file')
  info[1] <- paste('best gamma(s)',paste(bestgammas,collapse='\t'))
  info[2] <- paste('median of fl',paste(fl.allModels[[ctr]],collapse="\t"))
  info[3] <- paste('median of nl',paste(nl.allModels[[ctr]],collapse="\t"))
  writeLines(text=as.character(info),paste(outPath,'/details.txt',sep=''))
  print('done writing the details file')
  
  ##########################################################################################
  #####################      plot coef heatmap      ########################################
  ##########################################################################################
  pdf(paste(outPath,'fusedLasso_coefs_heatmap',fileName,'.pdf',sep=''))
  print(fl$cv.fl$bestsol$beta)
  plot.histone.coef(apply(cv.beta.mat.fl,2,FUN=median),bin.cnt,feature.names,main=paste(fileName,'median of trials',sep='_'))#,cluster_rows = F, cluster_cols = F)
  plot.histone.coef(apply(cv.beta.mat.fl,2,FUN=median),bin.cnt,feature.names,main=paste(fileName,'mean of trials',sep='_'))#,cluster_rows = F, cluster_cols = F)
  plot.histone.coef(fl$cv.fl$bestsol$beta,bin.cnt,feature.names,main=paste(fileName,'one model','gamma*',bestgammas[bestGammasIdx],sep='_'))#,cluster_rows = F, cluster_cols = F)
  dev.off()

  ##########################################################################################
  #####################      scatter plots test data      ##################################
  ##########################################################################################
  cv.beta.mat.fl <- fl$cv.fl$bestsol$cv.beta.mat
  pdf(paste(outPath,'scatterPlots_test_',fileName,'.pdf',sep=''))
  plot.scatter(group.rescale(partition$test$x,0,1,bin.cnt,min(min(partition$train$x)),max(max(partition$train$x))),partition$test$y,fl$cv.fl$bestsol$beta,0,ylab,main=paste('fl',fileName,sep='_'),xlab='prediction')
  plot.scatter(group.rescale(partition$test$x,0,1,bin.cnt,min(min(partition$train$x)),max(max(partition$train$x))),partition$test$y,apply(cv.beta.mat.fl,2,FUN=median),0,ylab,main=paste('fl',fileName,'median',sep='_'),xlab='prediction')
  plot.scatter(group.rescale(partition$test$x,0,1,bin.cnt,min(min(partition$train$x)),max(max(partition$train$x))),partition$test$y,apply(cv.beta.mat.fl,2,FUN=mean),0,ylab,main=paste('fl',fileName,'mean',sep='_'),xlab='prediction')
#  plot.scatter(group.rescale(partition$test$x,0,1,bin.cnt,min(min(partition$train$x)),max(max(partition$train$x))),partition$test$y,fl$cv.fl,0.2,ylab,main=paste('fl',fileName,'outliersRemoved',sep='_'),xlab='prediction')
  ####################################################################

  plot.scatter(partition$test$x,partition$test$y,cv.beta.mat.nl$best.nl,0,ylab,main=paste('nl',fileName,sep='_'),xlab='prediction')
  plot.scatter(partition$test$x,partition$test$y,apply(cv.beta.mat.nl$cv.beta,2,FUN=median),0,ylab,main=paste('nl',fileName,'median',sep='_'),xlab='prediction')
  plot.scatter(partition$test$x,partition$test$y,apply(cv.beta.mat.nl$cv.beta,2,FUN=mean),0,ylab,main=paste('nl',fileName,'mean',sep='_'),xlab='prediction')
 # plot.scatter.nl(group.rescale(partition$test$x,0,1,bin.cnt,min(min(partition$train$x)),max(max(partition$train$x))),partition$test$y,nl.sh$cv.nl[[1]],0.2,ylab,main=paste('nl',fileName,'outliersRemoved',sep='_'),xlab='prediction')
  dev.off()
print('done test data scatter plots')
  ##########################################################################################
  #####################     scatter plots entire data      #################################
  ##########################################################################################
  pdf(paste(outPath,'scatterPlots_entireData_',fileName,'.pdf',sep=''))
  plot.scatter(group.rescale(rbind(partition$train$x,partition$test$x),0,1,bin.cnt,min(min(partition$train$x)),max(max(partition$train$x))),c(partition$train$y,partition$test$y),fl$cv.fl$bestsol$beta,0,ylab,main=paste('fl',fileName,sep='_'),xlab='prediction')
  plot.scatter(group.rescale(rbind(partition$train$x,partition$test$x),0,1,bin.cnt,min(min(partition$train$x)),max(max(partition$train$x))),c(partition$train$y,partition$test$y),apply(cv.beta.mat.fl,2,FUN=median),0,ylab,main=paste('fl',fileName,'median',sep='_'),xlab='prediction')
  plot.scatter(group.rescale(rbind(partition$train$x,partition$test$x),0,1,bin.cnt,min(min(partition$train$x)),max(max(partition$train$x))),c(partition$train$y,partition$test$y),apply(cv.beta.mat.fl,2,FUN=mean),0,ylab,main=paste('fl',fileName,'mean',sep='_'),xlab='prediction')

####################################################################

  plot.scatter(rbind(partition$train$x,partition$test$x),c(partition$train$y,partition$test$y),cv.beta.mat.nl$best.nl,0,ylab,main=paste('nl',fileName,sep='_'),xlab='prediction')
  plot.scatter(rbind(partition$train$x,partition$test$x),c(partition$train$y,partition$test$y),apply(cv.beta.mat.nl$cv.beta,2,FUN=median),0,ylab,main=paste('nl',fileName,'median',sep='_'),xlab='prediction')
  plot.scatter(rbind(partition$train$x,partition$test$x),c(partition$train$y,partition$test$y),apply(cv.beta.mat.nl$cv.beta,2,FUN=mean),0,ylab,main=paste('nl',fileName,'mean',sep='_'),xlab='prediction')
  dev.off()
  print('done entire data scatter plots')
  ##########################################################################################
  ##########################################################################################
  ##########################################################################################
  beta.mat.med <- matrix(apply(cv.beta.mat.fl,2,FUN=median)[2:ncol(cv.beta.mat.fl)],nrow=graphComponents,ncol=bin.cnt,byrow = T)
  beta.mat.mean <- matrix(apply(cv.beta.mat.fl,2,FUN=mean)[2:ncol(cv.beta.mat.fl)],nrow=graphComponents,ncol=bin.cnt,byrow = T)
  pdf(paste(outPath,'barplots_median_of_trials',fileName,'.pdf',sep=''))
  print(coefficients.barplot(beta.mat.med,'','sum over the bins',feature.names))
  print(bins.barplot(beta.mat.med,'Bins','sum over features',seq(-floor(bin.cnt/2)+1,floor(bin.cnt/2),1)))
  dev.off()

  pdf(paste(outPath,'barplots_mean_of_trials',fileName,'.pdf',sep=''))
  print(coefficients.barplot(beta.mat.mean,'','sum over the bins',feature.names))
  print(bins.barplot(beta.mat.mean,'Bins','sum over features',seq(-floor(bin.cnt/2)+1,ceiling(bin.cnt/2),1)))
  dev.off()
  ctr <- ctr + 1
  bestGammasIdx <- bestGammasIdx + 1
}
save(fl.allModels,nl.allModels,file=paste(outPath,'/fl_nl_allModels.RData',sep=''))
if(length(datasets)>1){
  pdf(paste(outPath,'/cross_models.pdf',sep=''))
  plot.fl.accuracy.allmodels(fl.allModels,datasets,datasets.names,main='flasso test set',percent=percent)
  plot.fl.accuracy.allmodels(nl.allModels,datasets,datasets.names,main='lasso test set',percent=percent)
  dev.off()
}
print('best gamma(s)')
print(bestgammas)
