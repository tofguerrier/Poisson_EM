#Validation EM
source("~/Dropbox/MscUCD_Christophe/R scripts/EM_CG_2014August.R")

test1 <- genMixPoi(n=200, MPI=c(0.8,0.2), 
                   lambda=matrix(c(1,50,1,50,2,60),nrow = 2, ncol = 3))
test1 <- genMixPoi(n=200, MPI=c(0.8,0.2), 
                   lambda=matrix(c(100,50,100,50),nrow = 2, ncol = 2))

data = test1$data
ztrue = test1$trueZ

bic <- EM_CG_Poisson(data = data,max_iterations = 5,max_group = 4, min_group = 4)

bic$pi
bic$membership == ztrue

test1 <- genMixPoi(n=200, pi=c(0.5,0.3,0.2), 
                   lambda=matrix(c(3,2,1,10,15,5,1,10,20),nrow = 3, ncol = 3))

data = test1$data
ztrue = test1$trueZ

bic <- EM_CG_Poisson(data = data,max_iterations = 100)

bic$pi
bic$membership == ztrue
#How to handle mislabelling ???

setwd("~/Dropbox/MscUCD_Christophe/")
#May Data
load(file="./data/data_math10030.Rda")
load(file="./data/data_math10030A.Rda")

index_Avideos = grep("^V", colnames(math10030A))
index_lec = grep("^Total_Lecture_Attended", colnames(math10030A))
index_msc = grep("^SupportCount", colnames(math10030A))
index_tut = grep("^TotTutorial", colnames(math10030A))
index_vid = grep("^TotVideos", colnames(math10030A))

#index_counts = c(index_Avideos,index_lec, index_msc, index_tut )
index_counts = c(index_lec, index_msc, index_tut,index_vid )

#video.mclust <- Mclust(math10030A[,index_counts])

video.pclustA <- EM_CG_Poisson(data = as.matrix(math10030A[,index_Avideos]),max_iterations = 10,max_group = 20)
video.pclust4 <- EM_CG_Poisson(data = as.matrix(math10030A[,index_counts]),max_iterations = 10,max_group = 20)
video.pclust <- EM_CG_Poisson(data = as.matrix(math10030A[,c(index_lec, index_msc, index_tut)]),max_iterations = 10)

video.pclust$pi
video.pclust$membership

#As we increase the number of variables the nuumber of potential groups increase !
math10030A$membership <- video.pclust$membership

par(mfrow=c(2,2))
for ( i in 1:4) {
  barplot(table(math10030A[math10030A$membership==i,]$Grade),
                ylab="Frequency", main=paste0("Cluster ", i))
  }

par(mfrow=c(2,2))
for ( i in 1:4) {
  barplot(table(math10030A[math10030A$membership==i,]$Entry_score),
          ylab="Frequency", main=paste0("Cluster ", i))
}

par(mfrow=c(2,2))
for ( i in 1:4) {
  barplot(table(math10030A[math10030A$membership==i,]$Total_Lecture_Attended),
          ylab="Frequency", main=paste0("Cluster ", i))
}

