#Generate finite mixture poisson
genMixPoi <- function(n=200, MPI, lambda){
  G=dim(lambda)[1]
  K=dim(lambda)[2]
  
  data=matrix(0, n, K)
  #saving group membership
  ztrue=rep(0,n)
  for(i in 1:n)
  {
    ztrue[i]=c(1:G)%*%rmultinom(1,1, prob=MPI)
    for(k in 1:K)
    {
      data[i,k]=rpois(1,lambda[ztrue[i],k])
    }
  }
  return(list(data=as.matrix(data), trueZ = ztrue))
}

#E step
E_calZ <- function(x, MPI, lambda, classification=FALSE){
  #Inputs
  # x  : data ( matrix n observation x K variables )
  # MPI : mixing proportions suming to 1. vector of length G
  # lambda: means matrix G x K ( one for each variable in each of the G clusters)
  z <-matrix(0, nrow = nrow(x),ncol=length(MPI))
  for(obs in 1:nrow(x)){
    for(ig in 1:length(MPI)){
      z[obs,ig] <- MPI[ig] * prod(dpois(x= x[obs,],lambda =  lambda[ig,]))
#      if(z[obs,ig] == 0 ) {z[obs,ig] = 10e-46}
#      if(z[obs,ig] == 1 ) {z[obs,ig] = 1 - 10e-46}
    }
  }
  z <- apply(z, 2, function(x) x/rowSums(z))
  
  ## Smoothing prior to M-Step (0's set to 1e-10, 1's set to 1-1e-10)
  ## epsilon <- 1e-10
  ## maxcut <- 1-epsilon
  ## mincut <- epsilon
  ## z <- apply(z, 2, pmax, mincut)
  ## z <- apply(z, 2, pmin, maxcut);
  ## z <- z / rowSums(z)
  return(z)
}

#M step
M_calLambda <- function(x,z) {
  
  lambda<-matrix(0, nrow = ncol(z),ncol=ncol(x))
  for(g in 1:ncol(z)){
    for(k in 1:ncol(x)){
      #mu[g,k] <- sum(z[,g] * x[,k]) / nrow(x)
      lambda[g,k] <- sum(z[,g] * x[,k]) / sum(z[,g])
      #if(lambda[g,k]==0) { lambda[g,k] = 1e-10}
      if(lambda[g,k]==0) { warning("LAMBDA is ZERO !!",g,k)}
    }
  }
  return(lambda)
}

M_calPi <- function(z) {
  MPI <- colSums(z) / nrow(z)
  return(MPI)
}

# Initialisation through hclust
EM_init <- function(data, ng, method="hclust"){
  if( method != "hclust") { warning("Only hclust methode implemented.")}  
  #Make sure data is a matrix
  k = ncol(data)
  nobs = nrow(data)
  data = as.data.frame(data)
  
  #calculate distance
  #1 standardize
  data.scaled <- scale(data, center=FALSE, scale=TRUE)
  
  #2 Euclidean distance ( does it make sense for count )
  dist.data = dist(data.scaled, method="euclidean")
  
  #3 hierarchical clustering ward method
  hclust.data <- hclust(dist.data, method="ward")
  #plot(hclust.data)
  data$membership <- cutree(hclust.data, k=ng)
    
  #Get mixing proportions
  MPI.t <- table(data$membership) / length(data$membership)
  MPI <- as.vector(MPI.t)
  
  #Get mean for each member
  lambda <- matrix(0,nrow=ng,ncol=k)
  for(g in 1:ng){
    means.group <-  apply(subset(data, membership == g),2, mean)
    lambda[g,] <- as.vector(means.group)[1:k]
  }
  return(list(MPI=MPI,lambda=lambda))
  
}

#BIC
CalcBIC <- function(x=data,MPI=PI_,lambda=LA_,z=Z_){
  nobs = nrow(x)
  nva = dim(lambda)[2]
  ngrp = dim(lambda)[1]
  #number of variables
  nvar = length(MPI) - 1 + dim(lambda)[1]*dim(lambda)[2]
  loglik <- 0
  for(obs in 1:nobs){
    for(group in 1:ngrp){
      poilik <- dpois(x = as.vector(x[obs,]),lambda = lambda[group,], log = TRUE)
      #if(poilik > signif(.Machine$double.xmax, 6)) {warning("sigma-squared falls below threshold")}
      loglik_element <- z[obs,group] * (log(MPI[group]) +  sum(poilik))
      #if(is.infinite(loglik_element)) { 
      #  cat("overflow\n")
      #  cat("x:",as.vector(x[obs,]),"lambda",lambda[group,])
      #  cat("mpi",MPI[group],"\n")
      #  cat("poilik",poilik,"\n")
      #  #overflow when lambda goes to 0 !!!
      #}
      loglik <- loglik + loglik_element
    }
  }
  bic <- 2 * loglik - nvar * log(nobs) 
  return(list(loglik=loglik,bic=bic))
}

#EM Algorithm
EM_CG_Poisson <- function(data,max_iterations=50,thresold=1e-10,max_group=10,min_group=1){
  bic <- matrix(0,nrow=max_group,ncol=max_iterations) #to store BIC values
  for(NumberGroup in min_group:max_group) {
    init.values <- EM_init(data,ng=NumberGroup)
    #First Iteration
    Z_ <- E_calZ(x = data, MPI=init.values$MPI,lambda = init.values$lambda,classification = FALSE)
    LA_ <- M_calLambda(x = data, z = Z_)
    PI_ <- M_calPi(z = Z_)
    #Iterations
    for(it in 1:max_iterations){
      Z_ <- E_calZ(x = data, MPI=PI_, lambda= LA_, classification = FALSE)
      LA_ <- M_calLambda(x = data, z = Z_)
      PI_ <- M_calPi(z = Z_)
      #Calculate BIC?likelihood(need to be stored)
      bic[NumberGroup,it] <- CalcBIC(x=data,MPI=PI_,lambda=LA_,z=Z_)$loglik
      #Plot to see how membership change
      #plot(data,col=apply(Z_, 1, which.max))
      par(mfrow=c(2,2))
      #plot(PI_,main = paste("Pi Iterations ",it),ylab="PI")
      #plot(LA_[1,],main = paste("Lambda Iterations ",it),ylab="Lambda")
      #plot(LA_[2,],main = paste("Lambda Iterations ",it),ylab="Lambda")
      #plot(LA_[3,],main = paste("Lambda Iterations ",it),ylab="Lambda")
      #plot(LA_[4,],main = paste("Lambda Iterations ",it),ylab="Lambda")
      par(mfrow=c(1,1))
      #Check for monotonicity !!
      if(it>=2){
        cat("IT:",it,"NG",NumberGroup,"\n")
        cat("BIC",bic[NumberGroup,it],"\n")
        cat("pi",PI_,"\n")
        cat("la",LA_,"\n")
        #cat("z",Z_,"\n")
        if(bic[NumberGroup,it-1] - bic[NumberGroup,it] > 0) {
          ### warning("Issue with monotonicity on iteration:",it," for ",NumberGroup," clusters.")
        }
      }
      #cat("BIC:",bic[NumberGroup,it],"\n")
      #cat("Lambda:",LA_,"\n")
      #cat("MPI:",PI_," sum:",sum(PI_),"\n")
    }
  plot(y=bic[NumberGroup,],x=1:it,main=paste("BIC",NumberGroup,""),ylab="BIC",xlab="Iterations",type = "l")
  }
  bic.max <- apply(bic,1,max)
  plot(y=bic.max,x=1:max_group,main="BIC",ylab="BIC",xlab="Number of clusters",type = "l")
  best.bic = which.max(bic.max) #max bic indicate number of cluster giving the best fit.
  
  #Second pass to get the variable ( change that)
  NumberGroup = best.bic
  init.values <- EM_init(data,ng=NumberGroup)
  #First Iteration
  Z_1 <- E_calZ(x = data, MPI=init.values$MPI,lambda = init.values$lambda,classification = FALSE)
  LA_1 <- M_calLambda(x = data, z = Z_1)
  PI_1 <- M_calPi(z = Z_1)
  Z_ = Z_1
  LA_ = LA_1
  PI_ = PI_1
  #Iterations
  for(it in 1:max_iterations){
    Z_ <- E_calZ(x = data, MPI=PI_, lambda= LA_, classification = FALSE)
    LA_ <- M_calLambda(x = data, z = Z_)
    PI_ <- M_calPi(z = Z_)
  }
  membership = apply(Z_, 1, which.max)
  
  
  
  return(list(bic=bic.max,
              membership=membership,
              group=best.bic,
              lambda=LA_,MPI=PI_))
}


