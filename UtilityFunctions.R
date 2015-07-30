pkgTest <- function(pkg.name)
{
  if (!require(pkg.name,character.only = TRUE))
  {
    install.packages(x,dep=TRUE)
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

group.rescale <- function(x,a,b,bin.cnt){
  x.rescaled <- matrix(NA,nrow=nrow(x),ncol=ncol(x))
  groups.cnt <- ncol(x)/bin.cnt
  for(i in seq(groups.cnt)){
    group <- x[,seq((i-1)*bin.cnt+1,i*bin.cnt,1)]
    group.min <- min(min(group))
    group.max <- max(max(group))
    #x.rescaled[,seq((i-1)*bin.cnt+1,i*bin.cnt,1)] <- ((group-t(replicate(n=nrow(group),group.min)))*(b-a)/(group.max-group.min)+a)
    x.rescaled[,seq((i-1)*bin.cnt+1,i*bin.cnt,1)] <- ((group-(matrix(group.min,nrow=nrow(group),ncol=ncol(group))))*(b-a)/(group.max-group.min)+a)
    
  }
  return(x.rescaled)
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

get.rss <- function(pred,y){
  n <- length(y)
  return((sqrt(1/n*sum((pred-y)^2))))
}

cv.fusedlasso <- function(x,y,method=c("fusedlasso","fusedlasso1d","fusedlasso2d"),nfold=10,...){
  pkgTest("genlasso")
  pkgTest("foreach")
  n <- nrow(x)
  p <- ncol(x)
  foldsz <- floor(n/nfold)
  if(foldsz<2){
    print(paste('Number of fold for CV (',nfold,') is too large with respect to the training size... Setting it to 2',sep=''))
    nfold <- 2
    foldsz <- floor(n/nfold)
  }
  outidx <- foreach (i=seq(nfold)) %dopar% seq(((i-1)*foldsz + 1),i*foldsz,1)
  train <- foreach (i=seq(nfold)) %dopar% list(x=x[-outidx[[i]],],y=y[-outidx[[i]]])
  validation <- foreach (i=seq(nfold)) %dopar% list(x=as.matrix(x[outidx[[i]],]),y=y[outidx[[i]]])
  cv.fl <- foreach (i=seq(nfold),.packages = "genlasso") %dopar% do.call(method[1],list(X=train[[i]]$x,y=train[[i]]$y,...))
  best <- Inf
  ctr <- 1
  for(fl in cv.fl){
    predictions <- sapply(1:ncol(fl$fit),function(i) return(as.matrix(validation[[ctr]]$x) %*% fl$beta[,i]))
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
