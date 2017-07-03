source('UtilityFunctions_sequentialCV.R')
get.outerCV.partitions <-function(x,y,n.folds=10){
  fold.size <- ceiling(nrow(x)/n.folds)
  partitions <- list();
  for(i in seq(n.folds)){
    start <- (i-1)*fold.size+1;
    end <- i*fold.size;
    if(end > nrow(x))
      end <- nrow(x);
    test.idx <- seq(start,end)
    train.idx <- seq(nrow(x))[-test.idx]
    partitions[[i]] <- list(train=list(x=x[train.idx,],y=y[train.idx]),test=list(x=x[test.idx,],y=y[test.idx]))
  }
  return(partitions)
}
get.gammas <- function(gammas,best.gamma,fixedSteps=20){
		print(c('inside get.gammas',gammas))
		g.min <- min(gammas)
		g.max <- max(gammas)
		dist <- c((best.gamma - g.min),(g.max-best.gamma))
		best.dist <- min(dist)
		if(g.min == best.gamma){
				new.interval <- c(best.gamma,g.max/2)
		}else{
				new.interval <- c((best.gamma - g.min/2),(best.gamma + g.min/2))
				if((best.gamma - g.min/2) < 0)
						new.interval <- c(0,(best.gamma + g.min/2))
		}
		stepSize <- (new.interval[2] - new.interval[1])/fixedSteps
		print(c(new.interval,stepSize))
		if(new.interval[1] + stepSize > best.gamma)
				return(best.gamma)
		return(seq(new.interval[1],new.interval[2],stepSize))
}
fusedlasso.main <- function(x,y,bin.cnt,edgs,gammas){#in parallel
  ###To compute the fusedlasso for the given data
  ###Input:
  #####x: input space (the set of features)
  #####y: target (response) vector
  #####bin.cnt: the No of bins
  #####edgs: the edges to construct the adjacency graph for fusedlasso
  #####gammas: a numeric vector representing the gammas for fusedlasso
  #pkgTest("doSNOW")
  pkgTest("genlasso")
  pkgTest("parallel")
  pkgTest("scatterplot3d")

  clust.cnt <- nfold*length(gammas)
  print(clust.cnt)
  if(clust.cnt >= 30)
		  clust.cnt <- 30
  print(clust.cnt)
  
  #pkgTest("cluster")
  hist.cnt <- ceiling(ncol(x)/bin.cnt) #the No of histone modifications
  best.rss <- Inf #keeps the best ratio
  fl.best <- NULL #keeps the best fusedlasso model based on the ratio

  edgs <- edgs + 1
  edgs <- c(1,1,edgs)

#  if(is.null(graph))
  gr <- graph(edgs,directed = F)
  #x <- group.rescale(x,0,1,bin.cnt)
  prev.rss <- Inf;
  curr.rss <- 0;
  stop.thresh <- 10^(-3)
  total.gammas.len <- 0
  search.depth <- 1
  cl <- makeCluster(min(length(gammas),(detectCores()-1)))#,"FORK")
  while((prev.rss - curr.rss) > stop.thresh){
    clusterExport(cl,list("cv.fusedlasso","x","y","bin.cnt","nfold","gammas","maxsteps","my.minlam","pkgTest","group.rescale","gr"),envir=environment())
    fl.res <- parSapply(cl,gammas,function(gamma){return(cv.fusedlasso(x,y,bin.cnt,method="fusedlasso",nfold=nfold,gr,gamma=gamma,maxsteps = maxsteps,minlam=my.minlam))})
    for(gamma.ctr in seq(length(gammas))){
		  if(fl.res[,gamma.ctr]$bestsol$validationMSE < best.rss){
			best.rss <- fl.res[,gamma.ctr]$bestsol$validationMSE;
			gamma.best <- gammas[gamma.ctr];
            fl.best <- fl.res[,gamma.ctr]
		  }
    }
    best.gamma <- gamma.best
    prev.rss <- curr.rss
    curr.rss <- best.rss
    gammas <- get.gammas(gammas,best.gamma)
    print(paste('search.depth=',search.depth,sep=''))
    search.depth <- search.depth + 1
    if(length(gammas == 1) && gammas == best.gamma){
		  flag <- F
		  print(c('gammas=',gammas,'best.gamma',best.gamma))
		  print('single gamma set')
		  bestGamma <- best.gamma
		  #all.cols <- c(all.cols,'red')
		  break;
    }
  }
  stopCluster(cl)
  if(flag){
  	print('(prev.rss - curr.rss) < stop.thresh')
  }
  print(c("outPath:",outPath))
  save(fl.res,file=paste(outPath,'fl.res.RData',sep=''))
  print(class(fl.res))
  err.best <- best.rss
  if(flag)
	bestGamma <- gamma.best
  print(bestGamma)
  x <- group.rescale(x,0,1,bin.cnt)
  valid.lambda.range <- T
  new.maxsteps <- maxsteps
  while(valid.lambda.range){
	  fl <- do.call("fusedlasso",list(X=cbind(1,x),y=y,graph=gr,gamma=bestGamma,maxsteps = new.maxsteps,minlam=my.minlam))
  	  if(min(fl$lambda) < fl.best$bestsol$lambda && max(fl$lambda) > fl.best$bestsol$lambda){
			  beta.optimal <- coef.genlasso(fl,fl.best$bestsol$lambda)$beta
  			  inc.order.idx <- order(fl$lambda)
			  beta.inc.ordered <- fl$beta[,inc.order.idx]
			  valid.lambda.range <- F
  			  optimal.lambda.idx <- findInterval(fl.best$bestsol$lambda,fl$lambda[inc.order.idx])
	  }
	  new.maxsteps <- 2*new.maxsteps
  }
  print('in fusedloass.main')
  pred <- cbind(1,x) %*% beta.optimal
  validationRMSE <- get.rss(pred,y) 
  #bestsol <- list(lambda=fl.best$bestsol$lambda,beta=beta.inc.ordered[,optimal.lambda.idx],df=summary(fl)[optimal.lambda.idx,1],validationRMSE=validationRMSE,cv.beta.mat=fl.best$cv.beta.mat,cv.bestERR=err.best)
  bestsol <- list(lambda=fl.best$bestsol$lambda,beta=beta.optimal,df=summary(fl)[optimal.lambda.idx,1],validationRMSE=validationRMSE,cv.beta.mat=fl.best$cv.beta.mat,cv.bestERR=err.best)
  fl.best <- list(bestobj=fl,bestsol=bestsol)
  return(list(cv.fl=fl.best,gamma.best=bestGamma))
  for(gamma in gammas){
    print(paste('gamma: ',gamma,sep=''))
    print('Computing flasso...')
    cv.fl <- cv.fusedlasso(x,y,method="fusedlasso",nfold=nfold,gr,gamma=gamma,maxsteps = maxsteps)
    print('Done computing flasso')
    
    df <- ncol(x) - cv.fl$bestsol$df
    rss <- cv.fl$bestsol$validationMSE
    if(err.best > rss){
      err.best <- rss
      fl.best <- cv.fl
      bestGamma <- gamma
    }
  }
  x <- group.rescale(x,0,1,bin.cnt)
  fl <- do.call("fusedlasso",list(X=cbind(1,x),y=y,graph=gr,gamma=bestGamma,maxsteps = maxsteps))
  inc.order.idx <- order(fl$lambda)
  optimal.lambda.idx <- findInterval(fl.best$bestsol$lambda,fl$lambda[inc.order.idx])
  if(optimal.lambda.idx == 0){
    optimal.lambda.idx <- 1
  }
  print(c('optimal.lambda.idx',optimal.lambda.idx))
  beta.inc.ordered <- fl$beta[,inc.order.idx]
  print('in fusedloass.main')
  pred <- cbind(1,x) %*% beta.inc.ordered[,optimal.lambda.idx]
  validationRMSE <- get.rss(pred,y) 
  bestsol <- list(lambda=fl.best$bestsol$lambda,beta=beta.inc.ordered[,optimal.lambda.idx],df=summary(fl)[optimal.lambda.idx,1],validationRMSE=validationRMSE,cv.beta.mat=fl.best$cv.beta.mat)
  fl.best <- list(bestobj=fl,bestsol=bestsol)
  return(list(cv.fl=fl.best,gamma.best=bestGamma))
}
fusedlasso.main_ForLambdaInterpolation <- function(x,y,bin.cnt,edgs,gammas,intercept_mode){#in parallel
  ###To compute the fusedlasso for the given data
  ###Input:
  #####x: input space (the set of features)
  #####y: target (response) vector
  #####bin.cnt: the No of bins
  #####edgs: the edges to construct the adjacency graph for fusedlasso
  #####gammas: a numeric vector representing the gammas for fusedlasso
  pkgTest("doSNOW")
  pkgTest("scatterplot3d")
  pkgTest("parallel")
  pkgTest("genlasso")

  clust.cnt <- nfold*length(gammas)
  print(clust.cnt)
  if(clust.cnt >= 30)
	  clust.cnt <- 30
  print(clust.cnt)
  hist.cnt <- ceiling(ncol(x)/bin.cnt) #the No of histone modifications
  err.best <- Inf
  fl.best <- NULL #keeps the best fusedlasso model

  if(intercept_mode){
    edgs <- edgs + 1
    edgs <- c(1,1,edgs)
  }


  gr <- graph(edgs,directed = F)
  all.lens <- NULL; all.lambdas <- NULL
  all.rss <- NULL
  all.cols <- NULL
  best.rss <- Inf
  gamma.ctr <- 1
  flag <- T
  all.cvms <- NULL; all.gammas.info <- NULL;all.gammas <- NULL;
  print(length(gammas))
#  clusterExport(cl=cluster,varlist=c('all.lens','all.lambdas','all.cvms','all.gammas.info','all.gammas','all.rss','best.rss'))
  #fl.res <- parLapply(cluster,gammas,function(gamma){return(cv.fusedlasso(x,y,method="fusedlasso",nfold=nfold,gr,gamma=gamma,maxsteps = maxsteps))})
  total.gammas.len <- 0
  search.depth <- 1
  cluster <- makeCluster(min(length(gammas),max((detectCores() - 1),1)))#,"FORK")
  clusterExport(cluster,list("cv.beta.matrix.fl","cv.fusedlasso.interpolationOnLambdas","x","y","bin.cnt","nfold","gammas","maxsteps","my.minlam","pkgTest","group.rescale","gr"),envir=environment())
  fl.res <- parSapply(cluster,gammas,function(gamma){return(cv.fusedlasso.interpolationOnLambdas(x,y,bin.cnt,method="fusedlasso",nfold=nfold,gr,intercept_mode,gamma=gamma,maxsteps = maxsteps,minlam=my.minlam))})
  for(gamma.ctr in seq(length(gammas))){
		  lambda.table.ordered <- fl.res[,gamma.ctr]$lambda.table.ordered
		  cvm <- fl.res[,gamma.ctr]$cvm
		  all.gammas.info[[gamma.ctr]] <- lambda.table.ordered;
		  all.gammas <- c(all.gammas,rep(gammas[gamma.ctr],times=length(cvm)));
		  lambda.idx <- which.min(cvm);
		  if(cvm[lambda.idx] < best.rss){
				  best.rss <- cvm[lambda.idx];
				  gamma.best <- gammas[gamma.ctr];
				  lambda.best <- lambda.table.ordered[2,lambda.idx]
		  }
		  all.lens <- c(all.lens,length(cvm));
		  all.cvms <- c(all.cvms,cvm)
		  cols <- rep('black',times=length(cvm))
		  cols[which.min(cvm)] <- 'red'
		  all.cols <- c(all.cols,cols)
  }
  best.gamma <- gamma.best
  cv.fl <- fl.res[,which(gammas == best.gamma)]
  total.gammas.len <- total.gammas.len + length(gammas)
  gammas <- seq(best.gamma/2,3/2*best.gamma,by=.01)
  
  #gammas <- get.gammas(gammas,best.gamma)
  print(paste('search.depth=',search.depth,sep=''))
  search.depth <- search.depth + 1
  if(length(gammas == 1) && gammas == best.gamma)
  {
	flag <- F
	print(c('gammas=',gammas,'best.gamma',best.gamma))
	print('single gamma set')
	bestGamma <- best.gamma
    stopCluster(cluster)
    save(fl.res,file=paste(outPath,'fl.res.RData',sep=''))
    return(list(cv.fl=cv.fl,gamma.best=bestGamma,best.lambda=lambda.best))
  }
  ### If further search is required in gamma grid:
  fl.res <- parSapply(cluster,gammas,function(gamma){return(cv.fusedlasso.interpolationOnLambdas(x,y,bin.cnt,method="fusedlasso",nfold=nfold,gr,intercept_mode,gamma=gamma,maxsteps = maxsteps,minlam=my.minlam))})
  for(gamma.ctr in seq(length(gammas))){
		  lambda.table.ordered <- fl.res[,gamma.ctr]$lambda.table.ordered
		  cvm <- fl.res[,gamma.ctr]$cvm
		  all.gammas.info[[gamma.ctr]] <- lambda.table.ordered;
		  all.gammas <- c(all.gammas,rep(gammas[gamma.ctr],times=length(cvm)));
		  lambda.idx <- which.min(cvm);
		  if(cvm[lambda.idx] < best.rss){
				  best.rss <- cvm[lambda.idx];
				  gamma.best <- gammas[gamma.ctr];
				  lambda.best <- lambda.table.ordered[2,lambda.idx]
		  }
		  all.lens <- c(all.lens,length(cvm));
		  all.cvms <- c(all.cvms,cvm)
		  cols <- rep('black',times=length(cvm))
		  cols[which.min(cvm)] <- 'red'
		  all.cols <- c(all.cols,cols)
  }
  bestGamma <- gamma.best
  cv.fl <- fl.res[,which(gammas == best.gamma)]
  if(F){
  for(i in seq(length(all.gammas.info)))all.lambdas <- c(all.lambdas,all.gammas.info[[i]][2,])
  for(i in seq(length(all.gammas.info)))all.rss <- c(all.rss,all.gammas.info[[i]][3,])
  print(paste("lambda.best=",lambda.best))
  save(gammas,total.gammas.len,search.depth,all.lambdas,all.gammas,all.cvms,all.cols,all.gammas.info,file=paste(outPath,'optimizationPath.RData',sep=''))
  png(paste(outPath,'regularizationPath_allGamma_1.png',sep=''))
   scatterplot3d(log10(all.lambdas),log10(all.gammas),log10(all.cvms),pch=20,color=all.cols,box=F)
  dev.off()
  png(paste(outPath,'regularizationPath_allGamma_2.png',sep=''))
   scatterplot3d(log10(all.gammas),log10(all.lambdas),log10(all.cvms),pch=20,color=all.cols,box=F)
  dev.off()
  png(paste(outPath,'regularizationPath_allGamma_3.png',sep=''))
   scatterplot3d(log10(all.cvms),log10(all.lambdas),log10(all.gammas),pch=20,color=all.cols,box=F)
  dev.off()
  png(paste(outPath,'regularizationPath_allGamma_4.png',sep=''))
   scatterplot3d(log10(all.lambdas),log10(all.cvms),log10(all.gammas),pch=20,color=all.cols,box=F)
  dev.off()
  png(paste(outPath,'regularizationPath_allGamma_5.png',sep=''))
   scatterplot3d(log10(all.gammas),log10(all.cvms),log10(all.lambdas),pch=20,color=all.cols,box=F)
  dev.off()
  png(paste(outPath,'regularizationPath_allGamma_6.png',sep=''))
   scatterplot3d(log10(all.cvms),log10(all.gammas),log10(all.lambdas),pch=20,color=all.cols,box=F)
  dev.off()
  }
  stopCluster(cluster)
  save(fl.res,file=paste(outPath,'fl.res.RData',sep=''))
  return(list(cv.fl=cv.fl,gamma.best=bestGamma,best.lambda=lambda.best))
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

	  cv.nl.lambdas <- mclapply(seq(trials),function(trial)cv.glmnet(x=group.rescale(partition[1,trial]$train$x,0,1,bin.cnt),y=partition[1,trial]$train$y,parallel=F)$lambda.min,mc.cores=20)
	  cv.nl <- mclapply(seq(trials),function(trial)glmnet(x=group.rescale(partition[1,trial]$train$x,0,1,bin.cnt),y=partition[1,trial]$train$y,lambda=cv.nl.lambdas[[trials]]),mc.cores=20)
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



plot.histone.coef <- function(beta,bin.cnt,feat.names,main.title,intercept_mode){
  pkgTest("gplots")
  
  if(class(beta)[1]=="list")
    if(!is.null(beta$beta))
	  beta <- beta$beta
    else stop("The beta variable must either be a list containing an element named beta or a numeric vector of beta coefficients of the regressioni model!")
  if(intercept_mode)
     beta <- beta[2:length(beta)] #excluding intercept
  x <- beta
  x.matrix <- t(sapply(seq(floor(length(x)/bin.cnt)),function(i) x[seq((i-1)*bin.cnt+1,i*bin.cnt,1)]))
    colnames(x.matrix) <- seq((-floor(bin.cnt)/2)+1,ceiling(bin.cnt/2),1)
    rownames(x.matrix) <- feat.names
	if(nrow(x.matrix) == 1)
		image(x.matrix)
	else
	  heatmap.2(x.matrix, dendrogram="none", Rowv=FALSE, Colv=FALSE,lwid = c(.5,.5,4.5,.3),lmat=rbind(c(4,4,3,0),c(0,2,1,0)),
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

plot.stability.var <- function(cv.beta.mat.fl,cv.beta.mat.nl,intercept_mode){
  if(intercept_mode)
    cv.beta.mat.fl <- cv.beta.mat.fl[,-1]#remove intercept
  cv.beta.mat.nl <- cv.beta.mat.nl[,-1]#remove intercept
  fl.var <- coef.variation(cv.beta.mat.fl)
  nl.var <- coef.variation(cv.beta.mat.nl)
  boxplot(list(FusedLASSO=fl.var,LASSO=nl.var),col=c('red','blue'),ylab='variance of non-zero CV coefs',main='Stability')
}

plot.stability <- function(cv.beta.mat.fl,cv.beta.mat.nl,bin.cnt,feature.name,intercept_mode){
  par(mfrow = c(3,2))
  if(intercept_mode)
    cv.beta.mat.fl <- cv.beta.mat.fl[,-1]#remove intercept
  cv.beta.mat.nl <- cv.beta.mat.nl[,-1]#remove intercept
  range.max <- c(max(min(cv.beta.mat.fl),min(cv.beta.mat.nl)),max(max(cv.beta.mat.fl),max(cv.beta.mat.nl)))
  sapply(seq(length(feature.name)),function(i){boxplot(cv.beta.mat.fl[,seq((i-1)*bin.cnt,i*bin.cnt,1)],
                                                   col='red',main=paste('fl',feature.name[i],sep='_'),ylab='coefficients',ylim=range.max)
	                                               boxplot(cv.beta.mat.nl[,seq((i-1)*bin.cnt,i*bin.cnt,1)],
                                                   col='blue',main=paste('nl',feature.name[i],sep='_'),ylim=range.max)
  })
  par(mfrow = c(1,1))
}

plot.stability_shuffling <- function(fl.shuffling.obj,nl.shuffling.obj,bin.cnt,histone.names){
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
  
  fl <- foreach(trial = seq(trials),.packages = 'genlasso',.export = c("cv.fusedlasso","pkgTest")) %dopar%  cv.fusedlasso(shuffle[[trial]]$x,shuffle[[trial]]$y,method="fusedlasso",edges=edgs,nfold=nfold,
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

plot.scatter <- function(x,y,cv.fl,outlier.thresh,ylab='measured expr.',is_fl,intercept_mode,...){
  if(!is_fl | intercept_mode)
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
  plot(pred[-outliers],y[-outliers],ylab=ylab,col='blue',...)
  reg <- lm(y~pred)
  #abline(reg,col='red')
  legend('topleft',legend = paste('s=',correlation.spearman,'_p=',correlation.pearson,'_RMSE=',round(RSS,2),sep=''))
  return(list(rss=RSS,cor=c(correlation.spearman,correlation.pearson)))
}

plot.scatter.nl <- function(x,y,cv.nl,outlier.thresh,ylab='measured expr.',...){
  pred <- predict(cv.nl,x)
  RSS <- get.rss(pred,y)
  outliers <- order((y-pred)^2,decreasing = T)[1:floor(outlier.thresh*length(y))]
  correlation.spearman <- round(cor(y[-outliers],pred[-outliers],method='spearman'),2)
  correlation.pearson <- round(cor(y[-outliers],pred[-outliers],method='pearson'),2)
  ####Plot the scatter plot
  plot(pred[-outliers],y[-outliers],ylab=ylab,col='blue',...)
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
  pred.fl <- foreach(i=seq(trials),.export = "predict.fl") %do% predict.fl(fl.shuffling.obj$cv.fl[[i]]$bestsol,x)
  RSS.fl <- foreach(i=seq(trials),.export = "get.rss",.combine = 'c') %do% get.rss(pred.fl[[i]],y)
  correlation.fl <- foreach(i=seq(trials),.combine = 'c') %do% cor(y,pred.fl[[i]],method='spearman')
  
  pred.nl <- foreach(i = seq(trials)) %do% predict(nl.shuffling.obj$cv.nl[[i]],x)
  RSS.nl <- foreach(i=seq(trials),.export = "get.rss",.combine = 'c') %do% get.rss(pred.nl[[i]],y)
  correlation.nl <- foreach(i=seq(trials),.combine = 'c') %do% cor(y,pred.nl[[i]],method='spearman')
  
  return(list(fl=list(rss=round(mean(RSS.fl),2),cor=round(mean(correlation.fl),2)),nl=list(rss=round(mean(RSS.nl),2),cor=round(mean(correlation.nl),2))))

}



plot.fl.accuracy.allmodels <- function(models,datasets,datasets.names,main=NA,percent=.8){
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
	      partition <- data.partition(x = datasets[[ds.idx]]$x,y = datasets[[ds.idx]]$y,percent)
  		  partition$test$x <- group.rescale(partition$test$x,0,1,bin.cnt,min(min(partition$train$x)),max(max(partition$train$x)))
          pred <- cbind(1,partition$test$x) %*% models[[model.idx]]
	      #print(c('pred',pred))
	      rss <- get.rss(pred,partition$test$y)
	      correlations[model.idx,ds.idx] <- cor(partition$test$y,pred,method='spearman')
          RSSs[model.idx,ds.idx] <- rss
	      }
    }
	    
  rownames(correlations) <- headers
  colnames(correlations) <- headers
  rownames(RSSs) <- headers
  colnames(RSSs) <- headers
  pheatmap(correlations,main=paste('Spearman correlation of trained models vs\n all combinations_',main,sep=''),cluster_rows = T, cluster_cols = T,display_numbers=F,fontsize_number=8,fontsize_row = 8,fontsize_col = 8)
  pheatmap(correlations,main=paste('Spearman correlation of trained models vs\n all combinations_',main,sep=''),cluster_rows = T, cluster_cols = T,display_numbers=T,fontsize_number=8,fontsize_row = 8,fontsize_col = 8)
  pheatmap(RSSs,main=paste('RMSE of trained models vs\n all combinations_',main,sep=''),cluster_rows = T, cluster_cols = T,display_numbers=T,fontsize_number=8,fontsize_row = 8,fontsize_col = 8)
			  
  return(list(rss=RSSs,cor=correlations))
}



plot.fl.accuracy.allmodels_shuffling <- function(models,datasets,datasets.names,shuffle.idx,main=NA,percent=.8){
  ###Input
  #####models: a list containing di_sense, di_antisense, con_sense,and con_antisense elements of cv.fl objects
  #####datasets: a list containing di_sense, di_antisense, con_sense,and con_antisense elements of datasets, dataset is itself  a list of x and y
  pkgTest("foreach")
  pkgTest("pheatmap")
  get.rss <- function(pred,y){
		  n <- length(y)
		  return((sqrt(1/n*sum((pred-y)^2))))
  }
  
  correlations <- matrix(ncol=length(datasets),nrow=length(datasets))
  RSSs <- matrix(ncol=length(datasets),nrow=length(datasets))
  headers <- datasets.names
  print('flasso: Computing cross-dataset performance...')
  for(model.idx in seq(length(datasets))){
    for(ds.idx in seq(length(datasets))){
      partition <- foreach(i=seq(length(models[[model.idx]]$cv.fl)),.export = "data.partition") %do%
        data.partition(x = datasets[[ds.idx]]$x[shuffle.idx[i,],],y = datasets[[ds.idx]]$y[shuffle.idx[i,]],percent)
      pred <- foreach(i=seq(length(models[[model.idx]]$cv.fl)),.combine = 'rbind',.export = c("predict.fl")) %do%
        t(predict.fl(models[[model.idx]]$cv.fl[[i]]$bestsol,cbind(1,partition[[i]]$test$x)))
      #print(c('pred',pred))
      rss <- mean(foreach(i=seq(length(models[[model.idx]]$cv.fl)),.combine = 'c',.export = c("get.rss")) %do%
                    get.rss(pred[i,],partition[[i]]$test$y))
      correlations[model.idx,ds.idx] <- mean(foreach(i=seq(length(models[[model.idx]]$cv.fl)),.combine = 'c') %do%
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
  pheatmap(correlations,main=paste('FL_Spearman correlation of trained models vs\n all combinations_',main,sep=''),cluster_rows = T, cluster_cols = T,display_numbers=T,fontsize_number=8,fontsize_row = 8,fontsize_col = 8)
  pheatmap(RSSs,main=paste('FL_RMSE of trained models vs\n all combinations_',main,sep=''),cluster_rows = T, cluster_cols = T,display_numbers=T,fontsize_number=8,fontsize_row = 8,fontsize_col = 8)
  
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
            t(predict.fl(models[[paste(grp.model,anchor.model,sep='_')]]$cv.fl[[i]]$bestsol,partition[[i]]$test$x))
          
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
  pheatmap(correlations,main=paste('FL_Spearman correlation of trained models vs\n all combinations_',main,sep=''),cluster_rows = T, cluster_cols = T,display_numbers=T,fontsize_number=8,fontsize_row = 8,fontsize_col = 8)
  pheatmap(RSSs,main=paste('FL_RSS of trained models vs\n all combinations_',main,sep=''),cluster_rows = T, cluster_cols = T,display_numbers=T,fontsize_number=8,fontsize_row = 8,fontsize_col = 8)
  
  return(list(rss=RSSs,cor=correlations))
}

plot.nl.accuracy.allmodels <- function(models,datasets,datasets.names,shuffle.idx,main=NA,percent=.8){
  ###Input
  #####models: a list containing di_sense, di_antisense, con_sense,and con_antisense elements of cv.fl objects
  #####datasets: a list containing di_sense, di_antisense, con_sense,and con_antisense elements of datasets, dataset is itself  a list of x and y
  pkgTest("foreach")
  pkgTest("pheatmap")
  pkgTest("glmnet")
  get.rss <- function(pred,y){
		  n <- length(y)
		  return((sqrt(1/n*sum((pred-y)^2))))
  }
  correlations <- matrix(ncol=length(datasets),nrow=length(datasets))
  RSSs <- matrix(ncol=length(datasets),nrow=length(datasets))
  headers <- datasets.names
  print('lasso: Computing cross-dataset performance...')
  print(c('length(models[[1]])',length(models[[1]])))
  for(model.idx in seq(length(datasets))){
    for(ds.idx in seq(length(datasets))){
      partition <- foreach(i=seq(length(models[[model.idx]]$cv.nl)),.export = "data.partition") %do%
        data.partition(x = datasets[[ds.idx]]$x[shuffle.idx[i,],],y = datasets[[ds.idx]]$y[shuffle.idx[i,]],percent)
      pred <- foreach(i=seq(length(models[[model.idx]]$cv.nl)),.combine = 'rbind') %do%
        t(predict(models[[model.idx]]$cv.nl[[i]],partition[[i]]$test$x))
      #print(c('pred',pred))
      rss <- mean(foreach(i=seq(length(models[[model.idx]]$cv.nl)),.combine = 'c',.export = c("get.rss")) %do%
                    get.rss(pred[i,],partition[[i]]$test$y))
      correlations[model.idx,ds.idx] <- mean(foreach(i=seq(length(models[[model.idx]]$cv.nl)),.combine = 'c',.export = c("get.rss")) %do%
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
pheatmap(correlations,main=paste('NL_Spearman correlation of trained models vs\n all combinations_',main,sep=''),cluster_rows = T, cluster_cols = T,display_numbers=T,fontsize_number=8,fontsize_row = 8,fontsize_col = 8)
pheatmap(RSSs,main=paste('NL_RMSE of trained models vs\n all combinations_',main,sep=''),cluster_rows = T, cluster_cols = T,display_numbers=T,fontsize_number=8,fontsize_row = 8,fontsize_col = 8)

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
library(grid)
bins.barplot <- function(beta.mat,x.name,y.name,feature.names){
		pkgTest("ggplot2")
		beta.df <- data.frame(feature.name=1:ncol(beta.mat),positive=sapply(seq(ncol(beta.mat)),function(i)sum(beta.mat[which(beta.mat[,i]>=0),i])),
								                        negative=sapply(seq(ncol(beta.mat)),function(i)sum(beta.mat[which(beta.mat[,i]<0),i])))
  
  beta.df.final <- data.frame(feature.name=rep(1:ncol(beta.mat),each=2),val=c(rbind(beta.df$positive,beta.df$negative)),status=rep(c('positive','negative'),times=ncol(beta.mat)))
    
    ggplot(beta.df.final) + 
	    aes(x = feature.name, y = val, fill = status) +
		    geom_bar(stat = "identity", position = "identity") + scale_x_continuous(labels=feature.names,name=x.name,breaks=1:length(feature.names)) +
			    scale_y_continuous(name=y.name) +
				    theme(axis.text.y = element_text(vjust = 0.5, hjust = 0)) + theme(plot.margin = unit(c(1,1,1,1), "cm")) + theme(axis.text.x = element_text(size=8,angle=90))
}
