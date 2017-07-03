pkgTest <- function(pkg.name)
{
  if (!require(pkg.name,character.only = TRUE))
  {
    install.packages(pkg.name,dep=TRUE)
    if(!require(pkg.name,character.only = TRUE))
      stop(paste("Package ",pkg.name," not found!"))
  }
}

get.fscore <- function(m1,m2){
  value <- ((m1$mse-m2$mse)/(m2$mse))/((m2$df-m1$df)/m2$df)
  #   p.value <- pf(-log(value),m1$df,m2$df)
  p.value <- 1-pf((value),m1$df,m2$df,log.p=F)
  return(list(value=value,p.value=p.value))
}

r.squared <- function(y,yhat){
  return(sum((yhat-mean(y))^2)/(sum((y-mean(y))^2)))
}

adj.r.squared <- function(r.squared,n,m){
  ### n: num of samples
  ### m: num of DF
  return(RsquareAdj(r.squared,n,m))
}

rescale <- function(x,a,b,featureTypes){
  library('foreach')
  x.min <- apply(x,2,FUN=min)
  x.max <- apply(x,2,FUN=max)
  stdr <- ((x-t(replicate(n=nrow(x),x.min)))*(b-a)/(x.max-x.min)+a)
  return(stdr)
}

group.rescale <- function(x,a,b,bin.cnt,group.min=NULL,group.max=NULL){
  x.rescaled <- matrix(NA,nrow=nrow(x),ncol=ncol(x))
  groups.cnt <- ncol(x)/bin.cnt
  if(is.null(group.min) && is.null(group.max)){
    for(i in seq(groups.cnt)){
      group <- x[,seq((i-1)*bin.cnt+1,i*bin.cnt,1)]
      group.min <- c(group.min,min(min(group)))
      group.max <- c(group.max,max(max(group)))
    }
  }
	
  for(i in seq(groups.cnt)){
    group <- x[,seq((i-1)*bin.cnt+1,i*bin.cnt,1)]
    #x.rescaled[,seq((i-1)*bin.cnt+1,i*bin.cnt,1)] <- ((group-t(replicate(n=nrow(group),group.min)))*(b-a)/(group.max-group.min)+a)
    x.rescaled[,seq((i-1)*bin.cnt+1,i*bin.cnt,1)] <- ((group-(matrix(group.min[i],nrow=nrow(group),ncol=ncol(group))))*(b-a)/(group.max[i]-group.min[i])+a)
  }
  return(list(x=x.rescaled,min=group.min,max=group.max))
}

rescale2 <- function(x,a,b,featureTypes){
  library('foreach')
  bin.cnt <- ncol(x)/featureTypes
  x.original <- x
  if(min(x)==max(x))
    return(0*x)
  res <- foreach(i = seq(featureTypes),.combine='cbind') %do%{
    x <- x.original[,seq((i-1)*bin.cnt+1,i*bin.cnt,1)]
    x.min <- apply(x,1,FUN=min)
    x.max <- apply(x,1,FUN=max)
    stdr <- ((x-x.min)*(b-a)/(x.max-x.min)+a)
  }
  return(res)
}

rescale2 <- function(x,a,b){
  x.min <- min(x)
  x.max <- max(x)
  if(x.min==x.max)
    return(0*x)
  
  return((x-x.min)*(b-a)/(x.max-x.min)+a)
}

data.shuffle <- function(x,y){
  #assumes that samples are along the rows, i.e. Nxp (N samples with p features)
  N <- dim(x)[1]
  shuffle.idx <- sample(N)
  x.shuffled <- x[shuffle.idx,]
  if(is.null(dim(y))){
    y.shuffled <- y[shuffle.idx]
  }else
    y.shuffled <- y[shuffle.idx,]
  return(list(x=x.shuffled,y=y.shuffled,idx=shuffle.idx))
}

data.partition <- function(x,y,percent){
  #partitions the whole dataset into percent% training set and (1-percent)% test set
  N <- dim(x)[1]
  train.cnt <- floor(N*percent)
  train <- list(x=x[1:train.cnt,],y=y[1:train.cnt])
  test <- list(x=x[(train.cnt+1):N,],y=y[(train.cnt+1):N])
  return(list(train=train,test=test))
}

predict.fl <- function(cv.obj,x){
  return(x %*% cv.obj$beta)
}

get.rss <- function(pred,y){ #computes the RMSE
  n <- length(y)
  return((sqrt(1/n*sum((pred-y)^2))))
}

coef.variation <- function(mat){
  zero.idxs <- which(colSums(mat)==0)
  if(length(zero.idxs) == 0)
		  return(apply(mat,2,FUN=var))
  return(apply(mat[,-zero.idxs],2,FUN=var))
}

cv.beta.matrix.nl <- function(x,y,nfold){
	pkgTest("glmnet")
	pkgTest("foreach")
	n <- nrow(x)
	p <- ncol(x)
	foldsz <- floor(n/nfold)
	if(foldsz<2){
	  print(paste('Number of fold for CV (',nfold,') is too large with respect to the training size... Setting it to 2',sep=''))
	  nfold <- 2
	  foldsz <- floor(n/nfold)
	}
	outidx <- foreach (i=seq(nfold)) %do% seq(((i-1)*foldsz + 1),i*foldsz,1)
	train <- foreach (i=seq(nfold)) %do% list(x=x[-outidx[[i]],],y=y[-outidx[[i]]])
	validation <- foreach (i=seq(nfold)) %do% list(x=as.matrix(x[outidx[[i]],]),y=y[outidx[[i]]])
	cv.gl.res <- cv.glmnet(x=x,y=y,nfold=nfold,keep=T)
	lambda.opt <- cv.gl.res$lambda.min
	beta.mat <- matrix(NA,nrow=nfold,ncol=(p+1))
	for(i in seq(nfold)){
	  fit <- glmnet(train[[i]]$x,train[[i]]$y,lambda=lambda.opt)
	  beta <- predict(fit,type='coef')
	  beta.mat[i,] <- as.numeric(beta)
	}
	nl <- glmnet(x,y,lambda=lambda.opt)
	nl.beta <- as.numeric(predict(nl,type='coef'))
	print(c('length(nl.beta)',length(nl.beta)))
	return(list(cv.beta=beta.mat,best.nl=nl.beta))
}


cv.beta.matrix.fl <- function(cv.fl,opt.lambda){
  beta.mat <- matrix(NA,nrow=length(cv.fl),ncol=nrow(cv.fl[[1]]$beta))
  i <- 1
  for(fl in cv.fl){
    inc.ord <- order(fl$lambda)
    hit <- findInterval(opt.lambda,fl$lambda[inc.ord])
    if(hit == 0)
      hit <- 1
      beta.mat[i,] <- fl$beta[,inc.ord[hit]]
      i <- i + 1
    }
  return(beta.mat)
}

stability.var.plot <- function(fl.var,nl.var){
  boxplot(list(FL=fl.var,NL=nl.var),col=c('red','blue'),ylab='variance of non-zero CV coefs',main='Stability')
}

beta.interpolate <- function(lambda,all.lambdas,beta){#all.lambdas must be sorted increasingly
  beta.new <- vector(mode='numeric',length=nrow(beta))
 beta.new[seq(1,nrow(beta),1)] <- sapply(seq(1,nrow(beta),1),function(i){		 
      hit <- findInterval(lambda,all.lambdas)
      if(hit <= 1)
        return(beta[1])
      else if((all.lambdas[hit] - all.lambdas[(hit - 1)]) == 0){
        slope <- beta[i,hit]
	  }else{
        slope <- (beta[i,hit] - beta[i,(hit - 1)]) / (all.lambdas[hit] - all.lambdas[(hit - 1)])
	  }
      return(slope * lambda) 
}
)
 if(1==2){
  for(i in seq(2,nrow(beta),1)){
      hit <- findInterval(lambda,all.lambdas)
      if(hit <= 1)
        beta.new[i] <- beta[1]
        else{
	        slope <- (beta[i,hit] - beta[i,(hit - 1)]) / (all.lambdas[hit] - all.lambdas[(hit - 1)])
	        beta.new[i] <- slope * lambda
	    }
    }
 }
    return(beta.new)
}

##used to be cv.fusedlasso.minOfall, but after discussing with Marcel, we decided to use this old implementation again
cv.fusedlasso <- function(x,y,bin.cnt,method=c("fusedlasso","fusedlasso1d","fusedlasso2d"),nfold=10,gr,...){
  pkgTest("genlasso")
  pkgTest("parallel")
  n <- nrow(x)
  p <- ncol(x)
  x <- group.rescale(x,0,1,bin.cnt)$x
  foldsz <- floor(n/nfold)
  if(foldsz<2){
    print(paste('Number of fold for CV (',nfold,') is too large with respect to the training size... Setting it to 2',sep=''))
    nfold <- 2
    foldsz <- floor(n/nfold)
  }
  outidx <- train <- validation <- cv.fl <- list()
  for(i in seq(nfold)) outidx[[i]] <- seq(((i-1)*foldsz + 1),i*foldsz,1)
  for(i in seq(nfold)) train[[i]] <- list(x=x[-outidx[[i]],],y=y[-outidx[[i]]])
  for(i in seq(nfold)) validation[[i]] <- list(x=as.matrix(x[outidx[[i]],]),y=y[outidx[[i]]])

  cl.cv <- makeCluster(mc <- getOption("cl.cores",length(train)))
  #parallel::clusterExport(cl=cl,varlist=c("fusedlasso","gr","method","genlasso","gammas","maxsteps","my.minlam"),envir = environment())
  parallel::clusterExport(cl=cl.cv,varlist=c("fusedlasso","gr","method","gammas","maxsteps","my.minlam"),envir = environment())
  cv.fl <- parLapply(cl.cv,train,function(tr)do.call(method[1],list(X=cbind(1,tr$x),y=tr$y,graph=gr,...)))
  stopCluster(cl.cv)
#  cv.fl <- foreach (i=seq(nfold),.packages = "genlasso") %dopar% do.call(method[1],list(X=train[[i]]$x,y=train[[i]]$y,graph=gr,...))
  print('finished running do.call(flasso,...)')
  ctr <- 1
  mses <- list()
  best.rss <- Inf
  cv.beta.mat <- NULL
  for(fl in cv.fl){
	lambda.cnt <- ncol(fl$fit)
    predictions <- sapply(1:lambda.cnt,function(i) return(as.matrix(cbind(1,validation[[ctr]]$x)) %*% fl$beta[,i]))
	mses[[ctr]] <- t(sapply(1:lambda.cnt,function(i) return(1/length(predictions[,i])*sum((predictions[,i]-validation[[ctr]]$y)^2))))
	cv.beta.mat <- rbind(cv.beta.mat,fl$beta[,which.min(mses[[ctr]])])
	if(min(mses[[ctr]]) < best.rss){
			best.rss <- min(mses[[ctr]])
			lambda.min.idx <- which.min(mses[[ctr]])
			best.fl <- fl
			optimal.beta <- fl$beta[,lambda.min.idx]
			validationMSE <- best.rss
	}
    ctr <- ctr + 1
  }
  df.optimal <- summary(best.fl)[lambda.min.idx,1]
  bestSol <- list(lambda=best.fl$lambda[lambda.min.idx],beta=cv.beta.mat,df=df.optimal,validationMSE=validationMSE)
  #save(lambda.table,lambda.table.ordered,cv.fl,cvm,all.beta.interpolated,all.cv.rss,file=paste(outPath,'/cv.flasso.test_minOfAll_',cv.fl[[1]]$call$gamma,'.RData',sep=''))
  print('done saving in cv.fusedlasso')
  return(list(bestobj=best.fl,bestsol=bestSol,cv.beta.mat=cv.beta.mat))
}
cv.fusedlasso.interpolationOnLambdas <- function(x,y,bin.cnt,method=c("fusedlasso","fusedlasso1d","fusedlasso2d"),nfold=10,gr,intercept_mode,...){
  pkgTest("genlasso")
  pkgTest("foreach")
  pkgTest("parallel")
  n <- nrow(x)
  p <- ncol(x)
#  if(intercept_mode)
#    x <- group.rescale(x,0,1,bin.cnt)$x
  foldsz <- floor(n/nfold)
  if(foldsz<2){
    print(paste('Number of fold for CV (',nfold,') is too large with respect to the training size... Setting it to 2',sep=''))
    nfold <- 2
    foldsz <- floor(n/nfold)
  }
  outidx <- train <- validation <- cv.fl <- list()
  for(i in seq(nfold)) outidx[[i]] <- seq(((i-1)*foldsz + 1),i*foldsz,1)
  for(i in seq(nfold)) train[[i]] <- list(x=x[-outidx[[i]],],y=y[-outidx[[i]]])
  for(i in seq(nfold)) validation[[i]] <- list(x=as.matrix(x[outidx[[i]],]),y=y[outidx[[i]]])

  #cl.cv <- parallel::makeCluster(mc <- getOption("cl.cores",length(train)))
  #parallel::clusterExport(cl=cl.cv,varlist=c("fusedlasso","gr","method","gammas","maxsteps","my.minlam"),envir = environment())
  if(intercept_mode){
    cv.fl <- lapply(train,function(tr)do.call(method[1],list(X=cbind(1,tr$x),y=tr$y,graph=gr,...)))
  }else{
    cv.fl <- lapply(train,function(tr)do.call(method[1],list(X=tr$x,y=tr$y,graph=gr,...)))
  }
  #stopCluster(cl.cv)
  #cv.fl <- lapply(seq(nfold),function(i)do.call(method[1],list(X=cbind(1,train[[i]]$x),y=train[[i]]$y,graph=gr,...)))
#  cv.fl <- foreach (i=seq(nfold),.packages = "genlasso") %dopar% do.call(method[1],list(X=train[[i]]$x,y=train[[i]]$y,graph=gr,...))
  print('finished running do.call(flasso,...)')
  ctr <- 1
  mses <- list()
  lambda.table <- NULL
  for(fl in cv.fl){
	  lambda.cnt <- ncol(fl$fit)
  if(intercept_mode){
    predictions <- sapply(1:lambda.cnt,function(i) return(as.matrix(cbind(1,validation[[ctr]]$x)) %*% fl$beta[,i]))
  }else{
    predictions <- sapply(1:lambda.cnt,function(i) return(as.matrix(validation[[ctr]]$x) %*% fl$beta[,i]))
  }
	  mses[[ctr]] <- t(sapply(1:lambda.cnt,function(i) return(1/length(predictions[,i])*sum((predictions[,i]-validation[[ctr]]$y)^2))))
	  lambda.table <- cbind(lambda.table,rbind(rep(ctr,times = lambda.cnt),fl$lambda,mses[[ctr]]))
    ctr <- ctr + 1
  }
  lambdas.cnt <- ncol(lambda.table)
  lambda.table.ordered <- lambda.table[,order(lambda.table[2,])]
  all.beta.interpolated <- NULL
  all.cv.rss <- matrix(NA,nrow=nfold,ncol=lambdas.cnt)
  fold.idx <- vector(length=nfold,mode='numeric')
  for(i in seq(lambdas.cnt)){
    all.cv.rss[lambda.table.ordered[1,i],i] <- lambda.table.ordered[3,i]
    rest.of.idxs <- seq(nfold)[-lambda.table.ordered[1,i]]
    for(j in rest.of.idxs){
	    hit.idx <- which(lambda.table.ordered[1,]==j)
      hit <- findInterval(lambda.table.ordered[2,i],lambda.table.ordered[2,hit.idx])[1];
      if(hit==0){ # the lambda is smaller than the range of lambdas searched in fold #j
  		  hit <- hit.idx[1]
  	   all.cv.rss[j,i] <- lambda.table.ordered[3,hit]
      }else{
#		  ord.inc <- order(cv.fl[[j]]$lambda)
		  #beta.interpolated <- beta.interpolate(lambda.table.ordered[2,i],cv.fl[[j]]$lambda[ord.inc],cv.fl[[j]]$beta[,ord.inc])
	  	print(class(cv.fl[[j]]))
		  print(class(lambda.table.ordered[2,i]))
	  	beta.interpolated <- coef.genlasso(cv.fl[[j]],lambda.table.ordered[2,i])$beta
		  all.beta.interpolated <- cbind(all.beta.interpolated,beta.interpolated)


#		  pred <- cv.fl[[j]]$call$X %*% beta.interpolated
          if(intercept_mode){
		    pred <- as.matrix(cbind(1,validation[[j]]$x)) %*% beta.interpolated
          }else{
		    pred <- as.matrix(validation[[j]]$x) %*% beta.interpolated
          }
		  #all.cv.rss[j,i] <- lambda.table.ordered[3,hit.idx[hit]]
		  all.cv.rss[j,i] <- 1/length(pred)*sum((pred-validation[[j]]$y)^2)
		  #all.cv.rss[j,i] <- get.rss(cv.fl[[j]]$call$y,pred)
      }
    }
	}
  cvm <- colMeans(all.cv.rss)
  #plot(log2(cvm+1),type='l',col='blue',main=cv.fl[[1]]$call$gamma)
  lambda.min.idx <- which.min(cvm)[1]
  best.model.idx <- which(cv.fl[[lambda.table.ordered[1,lambda.min.idx]]]$lambda==lambda.table.ordered[2,lambda.min.idx])[1]
  lambda.optimal <- lambda.table.ordered[2,lambda.min.idx]
  beta.optimal <- cv.fl[[lambda.table.ordered[1,lambda.min.idx]]]$beta[,best.model.idx]
  cv.beta.mat <- cv.beta.matrix.fl(cv.fl,lambda.optimal)
  df.optimal <- summary(cv.fl[[lambda.table.ordered[1,lambda.min.idx]]])[best.model.idx,1]
  bestSol <- list(lambda=lambda.optimal,beta=beta.optimal,df=df.optimal,validationMSE=cvm[lambda.min.idx])
  best.fl <- cv.fl[[lambda.table.ordered[1,lambda.min.idx]]]
  #save(lambda.table,lambda.table.ordered,cv.fl,cvm,all.beta.interpolated,all.cv.rss,file=paste(outPath,'/cv.flasso.test_',cv.fl[[1]]$call$gamma,'.RData',sep=''))
  print('done saving in cv.fusedlasso')
  return(list(bestobj=best.fl,bestsol=bestSol,cv.beta.mat=cv.beta.mat,lambda.table.ordered=lambda.table.ordered,cvm=cvm))
}

cv.fusedlasso_min_of_all_folds <- function(x,y,method=c("fusedlasso","fusedlasso1d","fusedlasso2d"),nfold=10,edges,...){
  pkgTest("genlasso")
  pkgTest("foreach")
  n <- nrow(x)
  p <- ncol(x)
  x <- group.rescale(x,0,1,bin.cnt)
  foldsz <- floor(n/nfold)
  if(foldsz<2){
    print(paste('Number of fold for CV (',nfold,') is too large with respect to the training size... Setting it to 2',sep=''))
    nfold <- 2
    foldsz <- floor(n/nfold)
  }
  outidx <- foreach (i=seq(nfold)) %dopar% seq(((i-1)*foldsz + 1),i*foldsz,1)
  train <- foreach (i=seq(nfold)) %dopar% list(x=x[-outidx[[i]],],y=y[-outidx[[i]]])
  validation <- foreach (i=seq(nfold)) %dopar% list(x=as.matrix(x[outidx[[i]],]),y=y[outidx[[i]]])


  edges <- edges + 1
  edges <- c(1,1,edges)


#  if(is.null(graph))
  gr <- graph(edges,directed = F)
  cv.fl <- foreach (i=seq(nfold),.packages = "genlasso") %dopar% do.call(method[1],list(X=cbind(1,train[[i]]$x),y=train[[i]]$y,graph=gr,...))
#  cv.fl <- foreach (i=seq(nfold),.packages = "genlasso") %dopar% do.call(method[1],list(X=train[[i]]$x,y=train[[i]]$y,graph=gr,...))
  print('finished running do.call(flasso,...)')
  best <- Inf
  ctr <- 1
  for(fl in cv.fl){
    predictions <- sapply(1:ncol(fl$fit),function(i) return(as.matrix(cbind(1,validation[[ctr]]$x)) %*% fl$beta[,i]))
    #predictions <- sapply(1:ncol(fl$fit),function(i) return(as.matrix(validation[[ctr]]$x) %*% fl$beta[,i]))
    #correlations <- sapply(1:ncol(fl$fit),function(i) return(cor(predictions[,i],validation[[ctr]]$y,'spearman')))
    mses <- sapply(1:ncol(fl$fit),function(i) return(1/length(predictions[,i])*sum((predictions[,i]-as.matrix(validation[[ctr]]$y))^2)))
    best.idx <- which.min(mses)
    if(mses[best.idx] < best){
      best.fl <- fl
      best <- mses[best.idx]
      bestSol <- list(lambda=fl$lambda[best.idx],beta=fl$beta[,best.idx],df=summary(fl)[best.idx,1],validationMSE=mses[best.idx])
    }
    ctr <- ctr + 1
  }
  return(list(bestobj=best.fl,bestsol=bestSol))
}
