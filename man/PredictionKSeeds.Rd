% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/VariousFunctions.R
\name{PredictionKSeeds}
\alias{PredictionKSeeds}
\title{PredictionKSeeds}
\usage{
PredictionKSeeds(test, Seed, num_clusters, A, numcolsideffects)
}
\arguments{
\item{test}{test drugs features matrix}

\item{Seed}{matrix of seeds initialize in the KSeed algorithm}

\item{num_clusters}{number of clusters desired}

\item{A}{matrix of Naive Bayes predictions scores, result of KSeedsScores function}

\item{numcolsideffects}{number of sideeffects}
}
\value{
predizioni matrix containing predictions for the various uncharacterized drugs
}
\description{
Function implementing predictions for uncharacterized drugs
}
\examples{
r <- 8
c <- 10
m0 <- matrix(0, r, c)
num_clusters=4
features<-apply(m0, c(1,2), function(x) sample(c(0,1),1))
#Generate a sample side effects binary matrix
r1 <- 8
c1 <- 15
m1 <- matrix(0, r1, c1)
side_effects<-apply(m1, c(1,2), function(x) sample(c(0,1),1))
folds<-CreateFolds(features,2)
i=0
train = features[folds != i,]
trainpharmat = side_effects[folds != i,]
test = features[folds == i,]
testpharmat = side_effects[folds == i,]
s<-RandomSeedGenerator(num_clusters,nrow(train))
Seed<-SeedSelection(train,num_clusters,s)
clusters<-KSeedsClusters (train,num_clusters,Seed,s)
A<-KSeedsScores(train,trainpharmat,num_clusters,Seed,s,clusters)
predizioni<-PredictionKSeeds(test,Seed,num_clusters,A,ncol(side_effects))
}

