#--------------------------------------------
#Samples theta matrix 

#Inputs:
# ntk J x K matrix containing the number of pixels allocated to each transect j and FST k
# ngroup maximum number of FSTs
# gamma1 parameter between 0 and 1 for the prior of the truncated stick-breaking prior
# ntransect number of transects in the data set

#Output:
# theta J x K matrix containing the containing the proportion of each FST k in each transect j

update.theta=function(ntk,ngroup,gamma1,ntransect){

    #sample v from a beta distribution
    n.greater.k=t(apply(ntk[,ngroup:1],1,cumsum))[,ngroup:1]
    
    #get theta from v1 using the stick-breaking equation 
    theta=v1=matrix(NA,ntransect,ngroup)
    tmp=rep(1,ntransect)
    for (i in 1:(ngroup-1)){
        v1[,i]=rbeta(ntransect,ntk[,i],n.greater.k[,i+1]+gamma1)
        theta[,i]=v1[,i]*tmp
        tmp=tmp*(1-v1[,i])
    }
    theta[,ngroup]=tmp    

    #output results
    theta
}

#--------------------------------------------
#Samples pi matrix

#Inputs:
# z P vector containing the FST currently assigned to each pixel p
# nmat P x H matrix containing the number of returns for each pixel p in each height class h
# N.minus.n P x H matrix containing the results of Nmat-nmat
# ngroup maximum number of FSTs
# nheight number of height classes

#Output:
# pi1 K x H matrix containing the containing the return probability of each FST k in each height class h

update.pi1=function(z,nmat,N.minus.n,ngroup,nheight){
  pi1=matrix(NA,ngroup,nheight)
  for (k in 1:ngroup){
    cond=z==k
    soma=sum(cond)
    if (soma>1){
      nmat1=colSums(nmat[z==k,])
      N.minus.n1=colSums(N.minus.n[z==k,])
    }  
    if (soma==1){
      nmat1=nmat[z==k,]
      N.minus.n1=N.minus.n[z==k,]
    }  
    if (soma==0){
      nmat1=N.minus.n1=rep(0,nheight)      
    } 
    pi1[k,]=rbeta(nheight,nmat1+1,N.minus.n1+1)
  }
  pi1
}

#--------------------------------------------
#Calculate log-likelihood

#Inputs:
# z P vector containing the FST currently assigned to each pixel p
# nmat P x H matrix containing the number of returns for each pixel p in each height class h
# Nmat P x H matrix containing the number of incoming light pulses for each pixel p in each height class h
# ngroup maximum number of FSTs
# nheight number of height classes
# pi1 K x H matrix containing the containing the return probability of each FST k in each height class h
# nloc number of pixels in the data

#Output:
#sum of log-likelihood

calc.llk=function(z,nmat,Nmat,ngroup,nheight,pi1,nloc){
  llk=matrix(NA,nloc,nheight)
  for (k in 1:ngroup){
    cond=z==k
    for (j in 1:nheight){
      llk[cond,j]=dbinom(nmat[cond,j],size=Nmat[cond,j],prob=pi1[k,j],log=T)  
    }
  }
  sum(llk)
}