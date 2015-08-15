source('UtilityFunctions.R')
maxsteps = 2000
fusedlasso.main <- function(x,y,bin.cnt,edgs,gammas){#in parallel
  ###To compute the fusedlasso for the given data
  ###Input:
  #####x: input space (the set of features)
  #####y: target (response) vector
  #####bin.cnt: the No of bins
  #####edgs: the edges to construct the adjacency graph for fusedlasso
  #####gammas: a numeric vector representing the gammas for fusedlasso
  pkgTest("doSNOW")
  pkgTest("genlasso")
  cluster = makeCluster(20)
  registerDoSNOW(cluster)
  
  pkgTest("cluster")
  hist.cnt <- ceiling(ncol(x)/bin.cnt) #the No of histone modifications
  err.df.ratio.best <- -Inf #keeps the best f-tests
  fl.best <- NULL #keeps the best fusedlasso model based on the ftest
#  gr <- graph(edgs,directed = F)
#  x <- scale(x)
  #x <- group.rescale(x,0,1,bin.cnt)
  for(gamma in gammas){
    print(paste('gamma: ',gamma,sep=''))
    print('Computing flasso...')
    cv.fl <- cv.fusedlasso(x,y,method="fusedlasso",nfold=10,edges=edgs,gamma=gamma,maxsteps = maxsteps)
    print('Done computing flasso')
    
    df <- ncol(x) - cv.fl$bestsol$df
    rss <- cv.fl$bestsol$validationMSE
    err.df.ratio <- df/rss
    if(err.df.ratio.best < err.df.ratio){
      err.df.ratio.best <- err.df.ratio
      fl.best <- cv.fl
      bestGamma <- gamma
    }
  }
  stopCluster(cluster)
  return(list(cv.fl=fl.best,gamma.best=bestGamma))
}

fusedlasso.shuffling.par <- function(x,y,fl,shuffle.idx,bin.cnt,trial=20,percent=.8,cor.ret=T,rss.ret=T){
  pkgTest("doSNOW")
  pkgTest("genlasso")
  #Shuffle the data "trials" times (in parallel)
  if(is.null(shuffle.idx))
    shuffle.idx <- foreach(trial = seq(trials),.combine = 'rbind') %dopar% sample(nrow(x))
  trials <- nrow(shuffle.idx)
  cluster = makeCluster(20)
  registerDoSNOW(cluster)

  gamma <- fl$gamma.best
  hist.cnt <- ncol(x)/bin.cnt
  partition <- foreach(trial = seq(trials),.export = "data.partition") %dopar% data.partition(x[shuffle.idx[trial,],],y[shuffle.idx[trial,]],percent)
  
  best.idx <- which(fl$cv.fl$bestobj$lambda==fl$cv.fl$bestsol$lambda)        
  cv.fl <- foreach(trial = seq(trials),.packages = 'genlasso',.export = c("cv.fusedlasso","pkgTest")) %dopar%  
    cv.fusedlasso(partition[[trial]]$train$x,partition[[trial]]$train$y,method="fusedlasso",edges=edgs,nfold=10,
                  gamma = gamma,maxsteps = fl$cv.fl$bestobj$call$maxsteps)
  if(cor.ret&rss.ret){
    pred.fl <- foreach(i=seq(trials),.export = "predict.fl") %dopar% predict.fl(cv.fl[[i]]$bestsol,partition[[i]]$test$x)
    RSS.fl <- foreach(i=seq(trials),.export = "get.rss",.combine = 'c') %dopar% get.rss(pred.fl[[i]],partition[[i]]$test$y)
    correlation.fl <- foreach(i=seq(trials),.combine = 'c') %dopar% cor(partition[[i]]$test$y,pred.fl[[i]],method='spearman')
    stopCluster(cluster)
    return(list(cv.fl=cv.fl,rss=mean(RSS.fl),cor=mean(correlation.fl)))
  }
  else if(cor.ret){
    pred.fl <- foreach(i=seq(trials),.export = "predict.fl") %dopar% predict.fl(cv.fl[[i]]$bestsol,partition[[i]]$test$x)
    correlation.fl <- foreach(i=seq(trials),.combine = 'c') %dopar% cor(partition[[i]]$test$y,pred.fl[[i]],method='spearman')
    stopCluster(cluster)
    return(list(cv.fl=cv.fl,cor=mean(correlation.fl)))
  }
  else if(rss.ret){
    pred.fl <- foreach(i=seq(trials),.export = "predict.fl") %dopar% predict.fl(cv.fl[[i]]$bestsol,partition[[i]]$test$x)
    RSS.fl <- foreach(i=seq(trials),.export = "get.rss",.combine = 'c') %dopar% get.rss(pred.fl[[i]],partition[[i]]$test$y)
    stopCluster(cluster)
    return(list(cv.fl=cv.fl,rss=mean(RSS.fl)))
  }
  stopCluster(cluster)
  return(list(cv.fl=cv.fl))
}


fusedlasso.shuffling <- function(x,y,fl,bin.cnt,shuffle.idx,trials=20,percent=.8,cor.ret=T,rss.ret=T){
  pkgTest("genlasso")
  print(bin.cnt)
  library(parallel)
  print('Done loading libraries.')
  #Shuffle the data "trials" times (in parallel)
  if(is.null(shuffle.idx))
    shuffle.idx <- t(sapply(seq(trials),function(i)sample(nrow(x))))
  trials <- nrow(shuffle.idx)
	
  gamma <- fl$gamma.best
  hist.cnt <- ncol(x)/bin.cnt
  partition <- sapply(seq(trials),function(trial)data.partition(as.matrix(x[shuffle.idx[trial,],]),y[shuffle.idx[trial,]],percent))
  best.idx <- which(fl$cv.fl$bestobj$lambda==fl$cv.fl$bestsol$lambda)
  print("Running lapply")
  
  cluster = makeCluster(20)
  registerDoSNOW(cluster)
  cv.fl <- lapply(seq(trials), function(trial){print(trial)
							      cv.fusedlasso(partition[1,trial]$train$x,partition[1,trial]$train$y,method="fusedlasso",edges=edgs,nfold=10,
												                  gamma = gamma,maxsteps = fl$cv.fl$bestobj$call$maxsteps)
})
  stopCluster(cluster)
  print("Done running lapply!")
  if(cor.ret&rss.ret){
    pred.fl <- sapply(seq(trials),function(i) predict.fl(cv.fl[[i]]$bestsol,cbind(1,group.rescale(partition[2,i]$test$x,0,1,bin.cnt,min(min(partition[1,i]$train$x)),max(max(partition[1,i]$train$x))))))
    RSS.fl <- sapply(seq(trials),function(i) get.rss(pred.fl[,i],partition[2,i]$test$y))
    correlation.fl <- sapply(seq(trials),function(i) cor(partition[2,i]$test$y,pred.fl[,i],method='spearman'))
    return(list(cv.fl=cv.fl,rss=mean(RSS.fl),cor=mean(correlation.fl)))
  }##The rest of the cases must also be corrected for the group.rescaling
  else if(cor.ret){
    pred.fl <- sapply(seq(trials),function(i) predict.fl(cv.fl[[i]]$bestsol,cbind(1,partition[2,i]$test$x)))
    correlation.fl <- sapply(seq(trials),function(i) cor(partition[2,i]$test$y,pred.fl[,i],method='spearman'))
    return(list(cv.fl=cv.fl,cor=mean(correlation.fl)))
  }
  else if(rss.ret){
    pred.fl <- sapply(seq(trials),function(i) predict.fl(cv.fl[[i]]$bestsol,cbind(1,partition[2,i]$test$x)))
    RSS.fl <- sapply(seq(trials),function(i) get.rss(pred.fl[,i],partition[2,i]$test$y))
    return(list(cv.fl=cv.fl,rss=mean(RSS.fl)))
  }
  return(list(cv.fl=cv.fl))
}

normalasso.shuffling <- function(x,y,bin.cnt,shuffle.idx,trials=20,percent=.8,cor.ret=T,rss.ret=T){
  pkgTest("glmnet")
  library(parallel)
  #Shuffle the data "trials" times (in parallel)
  if(is.null(shuffle.idx))
		      shuffle.idx <- t(sapply(seq(trials),function(i)sample(nrow(x))))
    trials <- nrow(shuffle.idx)
	  hist.cnt <- ncol(x)/bin.cnt
	  partition <- sapply(seq(trials),function(trial)data.partition(as.matrix(x[shuffle.idx[trial,],]),y[shuffle.idx[trial,]],percent))
	  print('after partiotioning')
	  cv.nl <- mclapply(seq(trials),function(trial)cv.glmnet(x=group.rescale(partition[1,trial]$train$x,0,1,bin.cnt),y=partition[1,trial]$train$y,parallel=F),mc.cores=20)
	  print('after mclapply')
	  if(cor.ret&rss.ret){
		    pred.nl <- sapply(seq(trials),function(i){x.rescaled <- group.rescale(partition[2,i]$test$x,0,1,bin.cnt,group.min=min(min(partition[1,i]$train$x)),group.max=max(max(partition[1,i]$train$x)));predict(cv.nl[[i]],x.rescaled)})
	  ##The rest of the cases must also be corrected for the group.rescaling
	  print('after pred.nl')
		    RSS.nl <- sapply(seq(trials), function(i) get.rss(pred.nl[,i],partition[2,i]$test$y))
		    correlation.nl <- sapply(seq(trials),function(i) round(cor(partition[2,i]$test$y,pred.nl[,i],method='spearman'),2))
		    return(list(cv.nl=cv.nl,rss=mean(RSS.nl),cor=mean(correlation.nl,na.rm=T)))
	  }
	  else if(cor.ret){
	      pred.nl <- sapply(seq(trials),function(i) predict(cv.nl[[i]],partition[2,i]$test$x))
	      correlation.nl <- sapply(seq(trials),function(i) round(cor(partition[2,i]$test$y,pred.nl[,i],method='spearman'),2))
	      return(list(cv.nl=cv.nl,cor=mean(correlation.nl)))
	  }
		  else if(rss.ret){
				      pred.nl <- sapply(seq(trials),function(i) predict(cv.nl[[i]],partition[2,i]$test$x))
		      RSS.nl <- sapply(seq(trials), function(i) get.rss(pred.nl[,i],partition[2,i]$test$y))
			      stopCluster(cluster)
			      return(list(cv.nl=cv.nl,rss=mean(RSS.nl)))
				    }
		    return(list(cv.nl=cv.nl))
}


normalasso.shuffling.par <- function(x,y,shuffle.idx,trial=20,percent=.8,cor.ret=T,rss.ret=T){
  pkgTest("doSNOW")
  pkgTest("glmnet")
  #Shuffle the data "trials" times (in parallel)
  if(is.null(shuffle.idx))
    shuffle.idx <- foreach(trial = seq(trials),.combine = 'rbind') %dopar% sample(nrow(x))
  trials <- nrow(shuffle.idx)

  cluster = makeCluster(20)
  registerDoSNOW(cluster)
  hist.cnt <- ncol(x)/bin.cnt
  partition <- foreach(trial = seq(trials),.export = "data.partition") %dopar% data.partition(as.matrix(x[shuffle.idx[trial,],]),y[shuffle.idx[trial,]],percent)
  
  cv.nl <- foreach(trial = seq(trials),.packages = "glmnet") %dopar% cv.glmnet(x=group.rescale(partition[[trial]]$train$x,0,1,bin.cnt),y=partition[[trial]]$train$y,parallel=T)
  las.coefs <- foreach(trial = seq(trials),.combine = 'rbind') %dopar% coef(cv.nl[[trial]])[2:length(coef(cv.nl[[trial]]))]
  if(cor.ret&rss.ret){
    pred.nl <- foreach(i = seq(trials)) %dopar% predict(cv.nl[[i]],group.rescale(partition[[i]]$test$x,0,1,bin.cnt,min(min(partition[[trial]]$train$x)),max(max(partition[[trial]]$train$x))))
    RSS.nl <- foreach(i=seq(trials),.export = "get.rss",.combine = 'c') %dopar% get.rss(pred.nl[[i]],partition[[i]]$test$y)
    correlation.nl <- foreach(i=seq(trials),.combine = 'c') %dopar% round(cor(partition[[i]]$test$y,pred.nl[[i]],method='spearman'),2)
    stopCluster(cluster)
    return(list(cv.nl=cv.nl,rss=mean(RSS.nl),cor=mean(correlation.nl)))
  }
  else if(cor.ret){
    pred.nl <- foreach(i = seq(trials)) %dopar% predict(cv.nl[[i]],partition[[i]]$test$x)
    correlation.nl <- foreach(i=seq(trials),.combine = 'c') %dopar% round(cor(partition[[i]]$test$y,pred.nl[[i]],method='spearman'),2)
    stopCluster(cluster)
    return(list(cv.nl=cv.nl,cor=mean(correlation.nl)))
  }
  else if(rss.ret){
    pred.nl <- foreach(i = seq(trials)) %dopar% predict(cv.nl[[i]],partition[[i]]$test$x)
    RSS.nl <- foreach(i=seq(trials),.export = "get.rss",.combine = 'c') %dopar% get.rss(pred.nl[[i]],partition[[i]]$test$y)
    stopCluster(cluster)
    return(list(cv.nl=cv.nl,rss=mean(RSS.nl)))
  }
  stopCluster(cluster)
  return(list(cv.nl=cv.nl))
}

plot.histone.coef <- function(beta,bin.cnt,feat.names,main.title){
  pkgTest("gplots")
  
  if(class(beta)[1]=="list")
		  if(!is.null(beta$beta))
				  beta <- beta$beta
		  else stop("The beta variable must either be a list containing an element named beta or a numeric vector of beta coefficients of the regressioni model!")
		  beta <- beta[2:length(beta)] #excluding intercept
  x <- beta
  x.matrix <- t(sapply(seq(floor(length(x)/bin.cnt)),function(i) x[seq((i-1)*bin.cnt+1,i*bin.cnt,1)]))
    colnames(x.matrix) <- seq((-floor(bin.cnt)/2)+1,ceiling(bin.cnt/2),1)
    rownames(x.matrix) <- feat.names
	  heatmap.2(x.matrix, dendrogram="none", Rowv=FALSE, Colv=FALSE,
				    col = bluered(256), scale="none", key=TRUE, density.info="none",
					    trace="none",  symm=F,symkey=T,symbreaks=T,colsep = 0:ncol(x.matrix),rowsep = 0:nrow(x.matrix),sepcolor = 'black',sepwidth = c(.001,.001),main=main.title)
	  

}
plot.histone.coef.pheatmap <- function(beta,bin.cnt,histone.name,...){
  ###To plot the coefficients obtained by fusedlasso for each histone modification
  ###Input
  #####beta: either a numeric vector or a list containing beta element to represent the coefficients values
  #####bin.cnt: the No of bins
  #####hist.name: a vector of histone names appeared in the order of their presence in the feature vector
  pkgTest("pheatmap")
  if(class(beta)[1]=="list")
    if(!is.null(beta$beta))
      beta <- beta$beta
  else stop("The beta variable must either be a list containing an element named beta or a numeric vector of beta coefficients of the regressioni model!")
  beta <- beta[2:length(beta)] #excluding intercept
  h <- length(histone.name)
  histones <- beta[seq(1,bin.cnt,1)]
  for(i in 2:h)
  {
    histones <- rbind(histones,beta[seq((i-1)*bin.cnt+1,(i)*bin.cnt,1)])
  }
  colNames <- seq(-bin.cnt/2 + 1,bin.cnt/2,1)
  rownames(histones) <- histone.name
  colnames(histones) <- colNames
  if(F){
    uniqueVals <- unique(beta)
    lens <- c(length(which(uniqueVals < 0)),length(which(uniqueVals > 0)))
    lens2 <- c(ceiling(length(which(uniqueVals < 0))/2),ceiling(length(which(uniqueVals > 0))/2))
    if(min(beta) < 0)
      bk = unique(c(seq(min(beta),min(beta)/lens[1],length=lens[1]),0,seq(max(beta)/lens[2],max(beta),length=lens[2])))
    else
      bk <- unique(c(0,seq(max(beta)/lens[2],max(beta),length=lens[2]-1)))
    col1 = colorRampPalette(c("navy", "white"))(lens[1]) #set the order of greys
    col3 = colorRampPalette(c("white", "firebrick3"))(lens[2]+1)
    colors2 <- unique(c(col1, col3))
    if(!min(bk))
      colors2 <- colorRampPalette(c("white", "firebrick3"))(length(bk))
  #   print(bk)
  #   print(colors2)
  #   print(beta)
  }
  pheatmap(histones,...)#,color=colors2,...)#,legend_breaks=bk,legend_labels=as.character(round(bk,3)),...)#quantile(bk,seq(0,1,length=length(bk)))
}

plot.stability <- function(fl.shuffling.obj,nl.shuffling.obj,bin.cnt,histone.names){
  trials <- min(length(fl.shuffling.obj$cv.fl),length(nl.shuffling.obj$cv.nl))
  hist.cnt <- length(histone.names)
  print(c('in plot.stability... trials=',trials))  
  las.coefs.fullBeta <- foreach(trial = seq(trials),.combine = 'rbind') %do% as.numeric(coef(nl.shuffling.obj$cv.nl[[trial]]))
  fl.coefs.fullBeta <- foreach(trial = seq(trials),.combine = 'rbind') %do% fl.shuffling.obj$cv.fl[[trial]]$bestsol$beta
  print(paste('in plot.stability',length(coef(nl.shuffling.obj$cv.nl[[trial]])) ))
  las.coefs <- foreach(trial = seq(trials),.combine = 'rbind') %do% as.numeric(coef(nl.shuffling.obj$cv.nl[[trial]]))[2:length(coef(nl.shuffling.obj$cv.nl[[trial]]))]
  fl.coefs <- foreach(trial = seq(trials),.combine = 'rbind') %do% fl.shuffling.obj$cv.fl[[trial]]$bestsol$beta[2:length(fl.shuffling.obj$cv.fl[[trial]]$bestsol$beta)]
  x <- seq(-(floor(bin.cnt)/2-1),floor(bin.cnt)/2,1)
  print(dim(las.coefs))
  print(dim(fl.coefs))
  par(mfrow=c(3,2))
  for(i in seq(hist.cnt))
  {
    histone.idx <- seq((i-1)*bin.cnt+1,i*bin.cnt,1)
    y.lim <- range(c(range(fl.coefs[,histone.idx]),range(las.coefs[,histone.idx])),na.rm = T)
    if(i < hist.cnt)
    {
      boxplot(fl.coefs[,histone.idx],names=x,main=paste('fl',histone.names[i],sep='_'),col='red',ylim = y.lim,ylab='coefficients')
      boxplot(las.coefs[,histone.idx],names=x,main=paste('nl',histone.names[i],sep='_'),col='blue',ylim = y.lim)
    }else{
      boxplot(fl.coefs[,histone.idx],names=x,main=paste('fl',histone.names[i],sep='_'),col='red',ylim = y.lim,ylab='coefficients',xlab='features')
      boxplot(las.coefs[,histone.idx],names=x,main=paste('nl',histone.names[i],sep='_'),col='blue',ylim = y.lim,xlab='features')
    }
  }
  return(list(fl=list(var=apply(fl.coefs.fullBeta,2,FUN=var),mean=apply(fl.coefs.fullBeta,2,FUN=mean),median=apply(fl.coefs.fullBeta,2,FUN=median)),
              nl=list(var=apply(las.coefs.fullBeta,2,FUN=var),mean=apply(las.coefs.fullBeta,2,FUN=mean),median=apply(las.coefs.fullBeta,2,FUN=median))))
}

plot.stability2 <- function(x,y,cv.fl,bin.cnt,trials,histone.names,shuffle.idx=NULL){
  ###To check the stability of coefficients on fusedlasso and normal lasso models
  ###Input:
  #####x,y: input and output variables for the model
  #####cv.fl: object returned by cv.fusedlasso function
  #####bin.cnt: No of bins
  #####trials: No of trials
  #####histone.names: a vector of histone names appeared in the order of their presence in the feature vector
  pkgTest("doSNOW")
  cluster = makeCluster(20)
  registerDoSNOW(cluster)
  pkgTest("glmnet")
  #Shuffle the data "trials" times (in parallel)
  if(is.null(shuffle.idx))
    shuffle.idx <- foreach(trial = seq(trials),.combine = 'rbind') %dopar% sample(nrow(x))
  gamma <- cv.fl$gamma.best
  trials <- nrow(shuffle.idx)
  hist.cnt <- ncol(x)/bin.cnt
  shuffle <- foreach(trial = seq(trials)) %dopar% list(x=x[shuffle.idx[trial,],],y=y[shuffle.idx[trial,]])
  train.glmnet <- foreach(trial = seq(trials),.packages = "glmnet") %dopar% cv.glmnet(x=shuffle[[trial]]$x,y=shuffle[[trial]]$y,parallel=T)
  las.coefs <- foreach(trial = seq(trials),.combine = 'rbind') %dopar% coef(train.glmnet[[trial]])[2:length(coef(train.glmnet[[trial]]))]
  best.idx <- which(cv.fl$cv.fl$bestobj$lambda==cv.fl$cv.fl$bestsol$lambda)        
  
  fl <- foreach(trial = seq(trials),.packages = 'genlasso',.export = c("cv.fusedlasso","pkgTest")) %dopar%  cv.fusedlasso(shuffle[[trial]]$x,shuffle[[trial]]$y,method="fusedlasso",edges=edgs,nfold=10,
                                                                                   gamma = gamma,maxsteps = cv.fl$cv.fl$bestobj$call$maxsteps)$bestobj

  fl.bestsol <- foreach(trial = seq(trials)) %dopar% list(lambda=fl[[trial]]$lambda[best.idx[1]],beta=fl[[trial]]$beta[,best.idx[1]])
  fl.coefs <- foreach(trial = seq(trials),.combine = 'rbind') %dopar% fl.bestsol[[trial]]$beta
  par(mfrow=c(3,2))
  for(i in seq(hist.cnt))
  {
    histone.idx <- seq((i-1)*bin.cnt+1,i*bin.cnt,1)
    y.lim <- range(c(range(fl.coefs[,histone.idx]),range(las.coefs[,histone.idx])),na.rm = T)
    if(i < hist.cnt)
    {
      boxplot(fl.coefs[,histone.idx],main=paste('fl',histone.names[i],sep='_'),col='red',ylim = y.lim,ylab='coefficients')
      boxplot(las.coefs[,histone.idx],main=paste('nl',histone.names[i],sep='_'),col='blue',ylim = y.lim)
    }else{
      boxplot(fl.coefs[,histone.idx],main=paste('fl',histone.names[i],sep='_'),col='red',ylim = y.lim,ylab='coefficients',xlab='features')
      boxplot(las.coefs[,histone.idx],main=paste('nl',histone.names[i],sep='_'),col='blue',ylim = y.lim,xlab='features')
    }
      
  }
  stopCluster(cluster)
  return(list(fl=list(var=apply(fl.coefs,2,FUN=var),mean=apply(fl.coefs,2,FUN=mean)),nl=list(var=apply(las.coefs,2,FUN=var),mean=apply(las.coefs,2,FUN=mean))))
}

###Prediction plots
#plot.scatter <- function(x,y,cv.fl,...){
#  pred <- predict.fl(cv.fl$bestsol,x)
#  RSS <- get.rss(pred,y)
#  correlation <- round(cor(y,pred,method='spearman'),2)
#  ####Plot the scatter plot
#  plot(pred,y,ylab='measured expr.',col='blue',...)
#  reg <- lm(y~pred)
#  abline(reg,col='red')
#  legend('topleft',legend = paste('cor=',correlation,sep=''))
#  return(list(rss=RSS,cor=correlation))
#}

plot.scatter <- function(x,y,cv.fl,outlier.thresh,...){
  x <- cbind(1,x)
  if(!is.list(cv.fl))
	  pred <- x %*% cv.fl
  else
	  pred <- predict.fl(cv.fl$bestsol,x)
  RSS <- get.rss(pred,y)

  outliers <- order((y-pred)^2,decreasing = T)[1:floor(outlier.thresh*length(y))]
  correlation.spearman <- round(cor(y[-outliers],pred[-outliers],method='spearman'),2)
  correlation.pearson <- round(cor(y[-outliers],pred[-outliers],method='pearson'),2)
  ####Plot the scatter plot
  plot(pred[-outliers],y[-outliers],ylab='measured expr.',col='blue',...)
  reg <- lm(y~pred)
  #abline(reg,col='red')
  legend('topleft',legend = paste('s=',correlation.spearman,'_p=',correlation.pearson,'_RMSE=',round(RSS,2),sep=''))
  return(list(rss=RSS,cor=c(correlation.spearman,correlation.pearson)))
}

plot.scatter.nl <- function(x,y,cv.nl,outlier.thresh,...){
  pred <- predict(cv.nl,x)
  RSS <- get.rss(pred,y)
  outliers <- order((y-pred)^2,decreasing = T)[1:floor(outlier.thresh*length(y))]
  correlation.spearman <- round(cor(y[-outliers],pred[-outliers],method='spearman'),2)
  correlation.pearson <- round(cor(y[-outliers],pred[-outliers],method='pearson'),2)
  ####Plot the scatter plot
  plot(pred[-outliers],y[-outliers],ylab='measured expr.',col='blue',...)
  reg <- lm(y~pred)
  #abline(reg,col='red')
  legend('topleft',legend = paste('s=',correlation.spearman,'_p=',correlation.pearson,'_RMSE=',round(RSS,2),sep=''))
  return(list(rss=RSS,cor=c(correlation.spearman,correlation.pearson)))
}
accuracy.comparison <- function(x,y,fl.shuffling.obj,nl.shuffling.obj){
  ###To compute the spearman correlation and RSS of fusedlasso and normal lasso models
#   
#   RSS.fl <- vector(mode='numeric',length = trials)
#   correlation.fl <- vector(mode='numeric',length = trials)
#   RSS.nl <- vector(mode='numeric',length = trials)
#   correlation.nl <- vector(mode='numeric',length = trials)
#   gr <- fl$cv.fl$bestobj$call$graph
#   mxstep <- fl$cv.fl$bestobj$call$maxsteps
  pkgTest("foreach")
  trials <- min(length(fl.shuffling.obj),length(nl.shuffling.obj))
  pred.fl <- foreach(i=seq(trials),.export = "predict.fl") %do% predict.fl(fl.shuffling.obj$cv.fl[[i]]$bestsol,cbind(1,x))
  RSS.fl <- foreach(i=seq(trials),.export = "get.rss",.combine = 'c') %do% get.rss(pred.fl[[i]],y)
  correlation.fl <- foreach(i=seq(trials),.combine = 'c') %do% cor(y,pred.fl[[i]],method='spearman')
  
  pred.nl <- foreach(i = seq(trials)) %do% predict(nl.shuffling.obj$cv.nl[[i]],x)
  RSS.nl <- foreach(i=seq(trials),.export = "get.rss",.combine = 'c') %do% get.rss(pred.nl[[i]],y)
  correlation.nl <- foreach(i=seq(trials),.combine = 'c') %do% cor(y,pred.nl[[i]],method='spearman')
  
  return(list(fl=list(rss=round(mean(RSS.fl),2),cor=round(mean(correlation.fl),2)),nl=list(rss=round(mean(RSS.nl),2),cor=round(mean(correlation.nl),2))))

}

plot.fl.accuracy.allmodels <- function(models,datasets,datasets.names,shuffle.idx,main=NA,percent=.8){
  ###Input
  #####models: a list containing di_sense, di_antisense, con_sense,and con_antisense elements of cv.fl objects
  #####datasets: a list containing di_sense, di_antisense, con_sense,and con_antisense elements of datasets, dataset is itself  a list of x and y
  pkgTest("foreach")
  pkgTest("pheatmap")
  correlations <- matrix(ncol=length(datasets),nrow=length(datasets))
  RSSs <- matrix(ncol=length(datasets),nrow=length(datasets))
  headers <- datasets.names
  print('flasso: Computing cross-dataset performance...')
  for(model.idx in seq(length(datasets))){
    for(ds.idx in seq(length(datasets))){
      partition <- foreach(i=seq(length(models[[model.idx]])),.export = "data.partition") %do%
        data.partition(x = datasets[[ds.idx]]$x[shuffle.idx[i,],],y = datasets[[ds.idx]]$y[shuffle.idx[i,]],percent)
      pred <- foreach(i=seq(length(models[[model.idx]])),.combine = 'rbind',.export = c("predict.fl")) %do%
        t(predict.fl(models[[model.idx]]$cv.fl[[i]]$bestsol,cbind(1,partition[[i]]$test$x)))
      #print(c('pred',pred))
      rss <- mean(foreach(i=seq(length(models[[model.idx]])),.combine = 'c',.export = c("get.rss")) %do%
                    get.rss(pred[i,],partition[[i]]$test$y))
      correlations[model.idx,ds.idx] <- mean(foreach(i=seq(length(models[[model.idx]])),.combine = 'c',.export = c("get.rss")) %do%
                                    cor(partition[[i]]$test$y,pred[i,],method='spearman'))
      #print(c(model.idx,ds.idx,foreach(i=seq(length(models[[model.idx]])),.combine = 'c',.export = c("get.rss")) %do%
                #cor(partition[[i]]$test$y,pred[i,],method='spearman')))
      RSSs[model.idx,ds.idx] <- rss
      #print(c(model.idx,ds.idx,correlations[model.idx,ds.idx]))
    }
  }
  #print(c('length(datasets)',length(datasets)))
  #print(c('headers',headers))

  rownames(correlations) <- headers
  colnames(correlations) <- headers
  rownames(RSSs) <- headers
  colnames(RSSs) <- headers
  #print(c('cor:',correlations))
  #print(c('RSSs:',RSSs))
  pheatmap(correlations,main=paste('FL_Spearman correlation of trained models vs\n all combinations_',main,sep=''),cluster_rows = F, cluster_cols = F,display_numbers=T,fontsize_number=18,fontsize_row = 18,fontsize_col = 18)
  pheatmap(RSSs,main=paste('FL_RSS of trained models vs\n all combinations_',main,sep=''),cluster_rows = F, cluster_cols = F,display_numbers=T,fontsize_number=18,fontsize_row = 18,fontsize_col = 18)
  
  return(list(rss=RSSs,cor=correlations))
}

plot.fl.accuracy.allmodels2 <- function(models,datasets,shuffle.idx,main=NA,percent=.8){
  ###Input
  #####models: a list containing di_sense, di_antisense, con_sense,and con_antisense elements of cv.fl objects
  #####datasets: a list containing di_sense, di_antisense, con_sense,and con_antisense elements of datasets, dataset is itself  a list of x and y
  pkgTest("foreach")
  pkgTest("pheatmap")
  correlations <- matrix(ncol=4,nrow=4)
  RSSs <- matrix(ncol=4,nrow=4)
  headers <- c('di_sense','di_antisense','con_sense','con_antisense')
  ii <- 1
  for(grp.model in c('di','con')){
    for(anchor.model in c('sense','antisense')){
      jj <- 1
      for(grp.dataset in c('di','con')){
        for(anchor.dataset in c('sense','antisense')){
          partition <- foreach(i=seq(length(models[[paste(grp.model,anchor.model,sep='_')]])),.export = "data.partition") %do%
            data.partition(x = datasets[[paste(grp.dataset,anchor.dataset,sep='_')]]$x[shuffle.idx[i,],],y = datasets[[paste(grp.dataset,anchor.dataset,sep='_')]]$y[shuffle.idx[i,]],percent)
          pred <- foreach(i=seq(length(models[[paste(grp.model,anchor.model,sep='_')]])),.combine = 'rbind',.export = c("predict.fl")) %do%
            t(predict.fl(models[[paste(grp.model,anchor.model,sep='_')]]$cv.fl[[i]]$bestsol,cbind(1,partition[[i]]$test$x)))
          
          rss <- mean(foreach(i=seq(length(models[[paste(grp.model,anchor.model,sep='_')]])),.combine = 'c',.export = c("get.rss")) %do%
            get.rss(pred[i,],partition[[i]]$test$y))
          correlations[ii,jj] <- mean(foreach(i=seq(length(models[[paste(grp.model,anchor.model,sep='_')]])),.combine = 'c',.export = c("get.rss")) %do%
            cor(partition[[i]]$test$y,pred[i,],method='spearman'))
          RSSs[ii,jj] <- rss
          jj <- jj + 1
        }
      }
      ii <- ii + 1
    }
  }
  rownames(correlations) <- headers
  colnames(correlations) <- headers
  rownames(RSSs) <- headers
  colnames(RSSs) <- headers
  pheatmap(correlations,main=paste('FL_Spearman correlation of trained models vs\n all combinations_',main,sep=''),cluster_rows = F, cluster_cols = F,display_numbers=T,fontsize_number=18,fontsize_row = 18,fontsize_col = 18)
  pheatmap(RSSs,main=paste('FL_RSS of trained models vs\n all combinations_',main,sep=''),cluster_rows = F, cluster_cols = F,display_numbers=T,fontsize_number=18,fontsize_row = 18,fontsize_col = 18)
  
  return(list(rss=RSSs,cor=correlations))
}

plot.nl.accuracy.allmodels <- function(models,datasets,datasets.names,shuffle.idx,main=NA,percent=.8){
  ###Input
  #####models: a list containing di_sense, di_antisense, con_sense,and con_antisense elements of cv.fl objects
  #####datasets: a list containing di_sense, di_antisense, con_sense,and con_antisense elements of datasets, dataset is itself  a list of x and y
  pkgTest("foreach")
  pkgTest("pheatmap")
  correlations <- matrix(ncol=length(datasets),nrow=length(datasets))
  RSSs <- matrix(ncol=length(datasets),nrow=length(datasets))
  headers <- datasets.names
  print('lasso: Computing cross-dataset performance...')
  for(model.idx in seq(length(datasets))){
    for(ds.idx in seq(length(datasets))){
      partition <- foreach(i=seq(length(models[[model.idx]])),.export = "data.partition") %do%
        data.partition(x = datasets[[ds.idx]]$x[shuffle.idx[i,],],y = datasets[[ds.idx]]$y[shuffle.idx[i,]],percent)
      pred <- foreach(i=seq(length(models[[model.idx]])),.combine = 'rbind',.export = c("predict.fl")) %do%
        t(predict(models[[model.idx]]$cv.nl[[i]],cbind(1,partition[[i]]$test$x)))
      #print(c('pred',pred))
      rss <- mean(foreach(i=seq(length(models[[model.idx]])),.combine = 'c',.export = c("get.rss")) %do%
                    get.rss(pred[i,],partition[[i]]$test$y))
      correlations[model.idx,ds.idx] <- mean(foreach(i=seq(length(models[[model.idx]])),.combine = 'c',.export = c("get.rss")) %do%
                                               cor(partition[[i]]$test$y,pred[i,],method='spearman'))
      RSSs[model.idx,ds.idx] <- rss
    }
  }
rownames(correlations) <- headers
colnames(correlations) <- headers
rownames(RSSs) <- headers
colnames(RSSs) <- headers
#print(c('cor:',correlations))
#print(c('RSSs:',RSSs))
pheatmap(correlations,main=paste('NL_Spearman correlation of trained models vs\n all combinations_',main,sep=''),cluster_rows = F, cluster_cols = F,display_numbers=T,fontsize_number=18,fontsize_row = 18,fontsize_col = 18)
pheatmap(RSSs,main=paste('NL_RSS of trained models vs\n all combinations_',main,sep=''),cluster_rows = F, cluster_cols = F,display_numbers=T,fontsize_number=18,fontsize_row = 18,fontsize_col = 18)

return(list(rss=RSSs,cor=correlations))
}

plot.nl.accuracy.allmodels2 <- function(models,datasets,shuffle.idx,main=NA,percent=.8){
  ###Input
  #####models: a list containing di_sense, di_antisense, con_sense,and con_antisense elements of cv.fl objects
  #####datasets: a list containing di_sense, di_antisense, con_sense,and con_antisense elements of datasets, dataset is itself  a list of x and y
  pkgTest("foreach")
  pkgTest("pheatmap")
  pkgTest("glmnet")
  correlations <- matrix(ncol=4,nrow=4)
  RSSs <- matrix(ncol=4,nrow=4)
  headers <- c('di_sense','di_antisense','con_sense','con_antisense')
  ii <- 1
  for(grp.model in c('di','con')){
    for(anchor.model in c('sense','antisense')){
      jj <- 1
      for(grp.dataset in c('di','con')){
        for(anchor.dataset in c('sense','antisense')){
          partition <- foreach(i=seq(length(models[[paste(grp.model,anchor.model,sep='_')]])),.export = "data.partition") %do%
            data.partition(x = datasets[[paste(grp.dataset,anchor.dataset,sep='_')]]$x[shuffle.idx[i,],],y = datasets[[paste(grp.dataset,anchor.dataset,sep='_')]]$y[shuffle.idx[i,]],percent)
          pred <- foreach(i=seq(length(models[[paste(grp.model,anchor.model,sep='_')]])),.combine = 'rbind') %do%
            t(predict(models[[paste(grp.model,anchor.model,sep='_')]]$cv.nl[[i]],partition[[i]]$test$x))
          
          rss <- mean(foreach(i=seq(length(models[[paste(grp.model,anchor.model,sep='_')]])),.combine = 'c',.export = c("get.rss")) %do%
                        get.rss(pred[i,],partition[[i]]$test$y))
          correlations[ii,jj] <- mean(foreach(i=seq(length(models[[paste(grp.model,anchor.model,sep='_')]])),.combine = 'c',.export = c("get.rss")) %do%
                                        cor(partition[[i]]$test$y,pred[i,],method='spearman'))
          RSSs[ii,jj] <- rss
          jj <- jj + 1
        }
      }
      ii <- ii + 1
    }
  }
  rownames(correlations) <- headers
  colnames(correlations) <- headers
  rownames(RSSs) <- headers
  colnames(RSSs) <- headers
  pheatmap(correlations,main=paste('NL_Spearman correlation of trained models vs\n all combinations_',main,sep=''),cluster_rows = F, cluster_cols = F,display_numbers=T,fontsize_number=18,fontsize_row = 18,fontsize_col = 18)
  pheatmap(RSSs,main=paste('NL_RSS of trained models vs\n all combinations_',main,sep=''),cluster_rows = F, cluster_cols = F,display_numbers=T,fontsize_number=18,fontsize_row = 18,fontsize_col = 18)
  
  return(list(rss=RSSs,cor=correlations))
}

coefficients.barplot <- function(beta.mat,x.name,y.name,feature.names){
		pkgTest("ggplot2")
		  beta.df <- data.frame(feature.name=1:nrow(beta.mat),positive=sapply(seq(nrow(beta.mat)),function(i)sum(beta.mat[i,which(beta.mat[i,]>=0)])),
								             negative=sapply(seq(nrow(beta.mat)),function(i)sum(beta.mat[i,which(beta.mat[i,]<0)])))
  
  positive.order <- order(beta.df$positive + beta.df$negative,decreasing = F)
    negative.order <- order(beta.df$negative,decreasing = F)
    
    beta.df.sorted <- beta.df[positive.order,]
	  
	  feature.names.sorted <- feature.names[positive.order]
	  beta.df.final <- data.frame(feature.name=rep(1:nrow(beta.mat),each=2),val=c(rbind(beta.df.sorted$positive,beta.df.sorted$negative)),status=rep(c('positive','negative'),times=nrow(beta.mat)))
	    
	    
	    ggplot(beta.df.final) + 
		    aes(x = feature.name, y = val, fill = status) +
			    geom_bar(stat = "identity", position = "identity") + scale_x_continuous(breaks=1:length(feature.names.sorted),labels=feature.names.sorted,name=x.name) +
				    scale_y_continuous(name=y.name,breaks=seq(floor(min(beta.df.final$val)),ceiling(max(beta.df.final$val)),0.5)) +
					    theme(axis.text.y = element_text(vjust = 0.5, hjust = 0)) +coord_flip()
}

bins.barplot <- function(beta.mat,x.name,y.name,feature.names){
		pkgTest("ggplot2")
		beta.df <- data.frame(feature.name=1:ncol(beta.mat),positive=sapply(seq(ncol(beta.mat)),function(i)sum(beta.mat[which(beta.mat[,i]>=0),i])),
								                        negative=sapply(seq(ncol(beta.mat)),function(i)sum(beta.mat[which(beta.mat[,i]<0),i])))
  
  beta.df.final <- data.frame(feature.name=rep(1:ncol(beta.mat),each=2),val=c(rbind(beta.df$positive,beta.df$negative)),status=rep(c('positive','negative'),times=ncol(beta.mat)))
    
    ggplot(beta.df.final) + 
	    aes(x = feature.name, y = val, fill = status) +
		    geom_bar(stat = "identity", position = "identity") + scale_x_continuous(labels=feature.names,name=x.name,breaks=1:length(feature.names)) +
			    scale_y_continuous(name=y.name) +
				    theme(axis.text.y = element_text(vjust = 0.5, hjust = 0))
}
