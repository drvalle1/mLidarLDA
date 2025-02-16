rm(list=ls())
library('readr')
library('Rcpp')
set.seed(140)

#read useful functions
setwd('/Volumes/Users/drvalle/GIT_models/mLidarLDA')
sourceCpp('code.cpp')
source('mLidarLDA_main.R')
source('mLidarLDA_functions.R')

#get simulated data
setwd('/Users/denisvalle/Documents/EBA stuff/simulated example')
nmat=data.matrix(read_csv('simdat nmat1.csv'))[,-1]
Nmat=data.matrix(read_csv('simdat NMAT.csv'))[,-1]
nloc=nrow(nmat);
nheight=ncol(nmat)

trans.info=as.data.frame(read_csv('simdat auxiliary.csv'))
ntransect=length(unique(trans.info$trans.id)); ntransect
max(trans.info$trans.id)
trans.id=trans.info$trans.id
  
#settings
gamma1=0.1
ngroup=5
ngibbs=1000
burnin=ngibbs*0.5

#run gibbs sampler
mod=mLidarLDA(nmat,Nmat,ngroup,ntransect,trans.id,
                  ngibbs,burnin,gamma1)

#look at convergence
plot(mod$logl[200:ngibbs],type='l')

#look at number of groups
boxplot(mod$MLE.theta)

#export results
setwd('/Users/denisvalle/Documents/EBA stuff/simulated example')
write_csv(as.data.frame(mod$logl),'res logl.csv')
write_csv(as.data.frame(mod$MLE.pi1),'res MLE pi1.csv')
write_csv(as.data.frame(mod$MLE.theta),'res MLE theta.csv')
