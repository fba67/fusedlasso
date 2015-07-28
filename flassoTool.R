args <- commandArgs(trailingOnly = T)
#options(warn = -1)
a <- length(args)
params <- readLines(args[1])
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
source('wrapper_functions.R')
datasets <- list()
fl.allModels <- list()
nl.allModels <- list()
trials <- 10 
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
  print(class(data.x))
  print(length(data.x))
  
  #data.x <- data.x
  data.y <- as.numeric(readLines(inPaths[(p+1)],))#[1:nrow(data.x)])
  print(data.y[1:20])
  data.y <- log2(as.numeric(readLines(inPaths[(p+1)],))+1)#[1:nrow(data.x)])
  print(data.y[1:20])
  print(class(data.y))
  print(paste('Done reading ',fileName,'...',sep=''))
  bin.cnt <- ncol(data.x)/graphComponents
  print(c(bin.cnt,graphComponents,ncol(data.x)))
  edgs <- NULL
  for(i in 1:graphComponents)
  {
    edgs <- c(edgs,bin.cnt*(i-1)+1, rep((bin.cnt*(i-1)+2):(bin.cnt*i - 1),each=2),bin.cnt*i)
  }
  datasets[[ctr]] <- list(x=data.x,y=data.y)
  partition <- data.partition(data.x,data.y,percent=.8)
  shuffle.idx <- matrix(nrow=trials,ncol=nrow(data.x))
  for(i in seq(trials))
    shuffle.idx[i,] <- sample(nrow(data.x))
  print(gamma)
  print(dim(partition$train$x))
  fl <- fusedlasso.main(partition$train$x,partition$train$y,bin.cnt,edgs,seq(gamma$start,gamma$end,gamma$step))
  
  bestgammas[bestGammasIdx] <- fl$gamma.best
  bestGammasIdx <- bestGammasIdx + 1
  if(all(fl$cv.fl$bestsol$beta==0)){
    print(bestgammas)
    stop('The trained Fussed lasso model contains all the coefficients as zero... The models is not trained properly (needs more parameter tuning, maybe!)')
  }
  print(bestgammas)
  print('Starting shuffling...')
  nl.sh <- normalasso.shuffling(data.x,data.y,shuffle.idx,trial=trials)#,cor.ret=T,rss.ret=T
  fl.sh <- fusedlasso.shuffling(data.x,data.y,fl,shuffle.idx,trial=trials,percent=.8)
  print('Done shuffling.')
  fl.allModels[[ctr]] <- fl.sh
  nl.allModels[[ctr]] <- nl.sh
  
  pdf(paste(outPath,'stability_',fileName,'.pdf',sep=''))
  stab <- plot.stability(fl.sh,nl.sh,bin.cnt,feature.names)
  dev.off()
  
  pdf(paste(outPath,'fusedLasso_coefs_heatmap',fileName,'.pdf',sep=''))
  print(fl$cv.fl$bestsol$beta)
  plot.histone.coef(stab$fl$median,bin.cnt,feature.names,main=paste(fileName,'median of trials',sep='_'),cluster_rows = F, cluster_cols = F)
  plot.histone.coef(stab$fl$mean,bin.cnt,feature.names,main=paste(fileName,'mean of trials',sep='_'),cluster_rows = F, cluster_cols = F)
  plot.histone.coef(fl$cv.fl$bestsol$beta,bin.cnt,feature.names,main=paste(fileName,'one model',sep='_'),cluster_rows = F, cluster_cols = F)
  dev.off()
  print(class(stab$fl$mean)) 
  save(fl,fl.sh,stab,file=paste(outPath,'/fl_shfld_',fileName,'.RData',sep=''))
  pdf(paste(outPath,'scaterPlots_test_',fileName,'.pdf',sep=''))
  plot.scatter(partition$test$x,partition$test$y,fl$cv.fl,0,main=paste('fl',fileName,sep='_'),xlab='prediction')
  plot.scatter(partition$test$x,partition$test$y,stab$fl$median,0,main=paste('fl',fileName,'median',sep='_'),xlab='prediction')
  plot.scatter(partition$test$x,partition$test$y,stab$fl$mean,0,main=paste('fl',fileName,'mean',sep='_'),xlab='prediction')
  plot.scatter(partition$test$x,partition$test$y,fl$cv.fl,0.2,main=paste('fl',fileName,'outliersRemoved',sep='_'),xlab='prediction')
  plot.scatter.nl(partition$test$x,partition$test$y,nl.sh$cv.nl[[1]],0,main=paste('nl',fileName,sep='_'),xlab='prediction')
  plot.scatter.nl(partition$test$x,partition$test$y,nl.sh$cv.nl[[1]],0.2,main=paste('nl',fileName,'outliersRemoved',sep='_'),xlab='prediction')
  dev.off()

  pdf(paste(outPath,'scaterPlots_entireData_',fileName,'.pdf',sep=''))
  plot.scatter(rbind(partition$train$x,partition$test$x),c(partition$train$y,partition$test$y),fl$cv.fl,0,main=paste('fl',fileName,sep='_'),xlab='prediction')
  plot.scatter.nl(cbind(partition$train$x,partition$test$x),c(partition$test$y,partition$test$y),nl.sh$cv.nl[[1]],0,main=paste('nl',fileName,sep='_'),xlab='prediction')
  dev.off()
  ctr <- ctr + 1
}
save(fl.allModels,nl.allModels,file=paste(outPath,'/fl_nl_allModels.RData',sep=''))
if(length(datasets)>1){
  pdf(paste(outPath,'/cross_models.pdf',sep=''))
  plot.fl.accuracy.allmodels(fl.allModels,datasets,datasets.names,shuffle.idx,main='flasso test set',percent=.8)
  plot.nl.accuracy.allmodels(nl.allModels,datasets,datasets.names,shuffle.idx,main='lasso test set',percent=.8)
  dev.off()
}
print('best gamma(s)')
print(bestgammas)
