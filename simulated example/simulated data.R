rm(list=ls())
library('ggplot2')
set.seed(1)

#settings
ntransects=100
npixels_trans=400
nheight=10
ngroup=2

#spatial coordinates of center of each transect
combo=expand.grid(x=1:10,y=1:10)
combo$trans.id=1:ntransects

#define theta matrix
combo$FST1=1
combo$FST2=0
cond=combo$y==5; combo[cond,c('FST1','FST2')]=cbind(rep(0,sum(cond)),
                                                    rep(1,sum(cond)))
cond=combo$y%in%c(4,6); combo[cond,c('FST1','FST2')]=c(0.5,0.5)
theta=combo[,c('FST1','FST2')]
rownames(theta)=paste0('t',1:ntransects)

#display theta matrix for FST1
ggplot()+
  geom_tile(data=combo,alpha=0.5,aes(x=x,y=y,fill=FST1))+
  scale_fill_gradient2(low = "cyan", mid = "red",high='purple',
                       limits=c(0,1),midpoint=0.5)

#display theta matrix for FST2
ggplot()+
  geom_tile(data=combo,alpha=0.5,aes(x=x,y=y,fill=FST2))+
  scale_fill_gradient2(low = "cyan", mid = "red",high='purple',
                       limits=c(0,1),midpoint=0.5)

#parameters
seq1=exp(-0.1*seq(from=1,to=10,length.out=nheight)); seq1[1]=0.6
seq2=exp(-0.4*seq(from=1,to=10,length.out=nheight)); seq2[1]=0.9
pi1=rbind(seq1,
          seq2)
rownames(pi1)=paste0('FST',1:2)
colnames(pi1)=paste0('h',1:nheight)

#show pi1 matrix
plot(1:nheight,pi1[1,],type='h',lwd=2,col='red',ylim=c(0,1))
for (i in 1:nheight){
  lines(rep(i+0.1,2),c(0,pi1[2,i]),col='blue',lwd=2)
}

#create combinations of transect id and pixel id
combo1=expand.grid(pixel.id=1:npixels_trans,
                  trans.id=1:ntransects)
ncombo=nrow(combo1)
combo1$order1=1:ncombo

#simulate data nmat
Nmat=matrix(rpois(ncombo*nheight,lambda=100),ncombo,nheight)
nmat=matrix(NA,ncombo,nheight)
for (i in 1:ncombo){
  #sample FST
  prob=theta[combo1$trans.id[i],]
  ind=rmultinom(1,size=1,prob)
  FST=which(ind==1)
  
  #sample nmat
  prob=pi1[FST,]
  nmat[i,]=rbinom(nheight,size=Nmat[i,],prob)
}

#export data
setwd('/Users/denisvalle/Documents/EBA stuff/simulated example')
NMAT=Nmat
colnames(nmat)=colnames(NMAT)=paste0('h',1:nheight)
rownames(nmat)=rownames(NMAT)=paste0('p',1:(ntransects*npixels_trans))
write.csv(nmat,'simdat nmat1.csv')
write.csv(NMAT,'simdat NMAT.csv')

#export auxiliary information
combo2=merge(combo,combo1,all=T); dim(combo1); dim(combo2)
combo3=combo2[order(combo2$order1),] #make sure order is correct
ind=which(colnames(combo3)%in%c('FST1','FST2','order1'))
write.csv(combo3[,-ind],'simdat auxiliary.csv',row.names=F)

#export true parameters
write.csv(theta,'simdat theta.csv')
write.csv(pi1,'simdat pi1.csv')
