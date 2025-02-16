rm(list=ls())
library('readr')
library('ggplot2')

#get auxiliary information
setwd('/Users/denisvalle/Documents/EBA stuff/simulated example')
trans.info=as.data.frame(read_csv('simdat auxiliary.csv'))
trans.info1=unique(trans.info[,c('x','y','trans.id')])

#get estimated values
pi1.estim=read.csv('res MLE pi1.csv')
theta.estim=read.csv('res MLE theta.csv')
colnames(theta.estim)=paste0('FST',1:ncol(theta.estim))
theta.estim1=cbind(trans.info1,theta.estim)

#get true parameter values
pi1.true=read.csv('simdata pi1.csv')[,-1]
theta.true=read.csv('simdata theta.csv')[,-1]

#display boxplot of estimated theta
boxplot(theta.estim)

#display estimated theta matrix for FST1
ggplot()+
  geom_tile(data=theta.estim1,alpha=0.5,aes(x=x,y=y,fill=FST1))+
  scale_fill_gradient2(low = "cyan", mid = "red",high='purple',
                       limits=c(0,1),midpoint=0.5)

#display estimated theta matrix for FST2
ggplot()+
  geom_tile(data=theta.estim1,alpha=0.5,aes(x=x,y=y,fill=FST2))+
  scale_fill_gradient2(low = "cyan", mid = "red",high='purple',
                       limits=c(0,1),midpoint=0.5)

#display estimated pi1
nheight=ncol(pi1.estim)
plot(1:nheight,pi1.estim[1,],type='h',col='red',ylim=c(0,1))
for (i in 1:nheight){
  lines(rep(i,2)+0.1,c(0,pi1.estim[2,i]),col='blue')
}
