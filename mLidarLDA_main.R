#Inputs:
# nmat P x H matrix containing the number of returns for each pixel p in each height class h
# Nmat P x H matrix containing the number of incoming light pulses for each pixel p in each height class h
# ngroup maximum number of FSTs
# ntransect number of transects in the data set
# trans.id transect identifier determining which pixel belong to which transect
# ngibbs number of iterations for the Gibbs sampler
# nburnin number of iterations to discard as burn-in
# gamma1 parameter between 0 and 1 for the prior of the truncated stick-breaking prior

#Outputs:
# pi1 Matrix containing posterior samples for the vertical height profile of each FST (rows are samples, columns are FST x height class combinations)
# theta Matrix containing posterior samples for the proportion of each FST in each transect (rows are samples, columns are FST x transect combinations)
# logl Vector that stores the log-likelihood (useful for determining model convergence)
# MLE.pi1 K x H matrix containing the MLE for return probability in each FST k in each height class h 
# MLE.theta J x K matrix containing the MLE for the proportion of each FST k in each transect j
# MLE.z P vector containing the MLE FST assigned to each pixel p
# MLE.iter vector with the iteration for which we obtained the MLE estimates
# last.z P vector containing the last results regarding the FST assigned to each pixel p

mLidarLDA=function(nmat,Nmat,ngroup,ntransect,trans.id,
                   ngibbs,burnin,gamma1){
    
    #useful pre-calculated quantities 
    nloc=nrow(nmat)
    nheight=ncol(nmat)
    N.minus.n=Nmat-nmat

    #initial parameter values
    z=sample(1:ngroup,size=nloc,replace=T)
    theta=matrix(1/ngroup,ntransect,ngroup)
    tmp=runif(ngroup*nheight)
    pi1=matrix(tmp,ngroup,nheight)

    #matrices to store results from gibbs sampler
    store.pi1=matrix(NA,ngibbs,nheight*ngroup)
    store.theta=matrix(NA,ngibbs,ngroup*ntransect)
    # store.z=matrix(NA,ngibbs,nloc) #too big
    store.logl=rep(NA,ngibbs)
    
    #useful settings
    norder=50 #re-order identified FSTs every norder iteration
    max.llk=-Inf #current maximum likelihood value

    #run gibbs sampler
    for (i in 1:ngibbs){
        print(i)
        
        #sample group allocation vector z
        ltheta=log(theta)
        lpi1=log(pi1)
        l1.minus.pi1=log(1-pi1)
        z=samplez(ltheta=ltheta,
                  nmat=nmat,
                  Nminusn=N.minus.n,
                  lpi1=lpi1,
                  l1minuspi1=l1.minus.pi1,
                  ngroup=ngroup,
                  nloc=nloc,
                  nheight=nheight,
                  randu=runif(nloc),
                  TransID=trans.id-1)

        #re-order groups if necessary
        if (i%%norder==0 & i<burnin){
            med=apply(theta,2,mean)
            order1=order(med,decreasing=T)
            theta=theta[,order1]
            pi1=pi1[order1,]
            znew=rep(NA,nloc)
            for (j in 1:ngroup){
                cond=z==order1[j]
                znew[cond]=j
            }
            z=znew
        }
        
        #sample pi1
        pi1=update.pi1(z=z,
                       nmat=nmat,
                       N.minus.n=N.minus.n,
                       ngroup=ngroup,
                       nheight=nheight)

        #calculate the number of locations in each transect 
        #assigned to each group
        ntk=calc_ntk(z=z-1,
                     ngroup=ngroup,nloc=nloc,ntransect=ntransect,
                     TransID=trans.id-1)
        
        #sample theta
        theta=update.theta(ntk=ntk,ngroup=ngroup,gamma1=gamma1,
                           ntransect=ntransect)

        #get loglikelihood
        logl=calc.llk(z=z,nmat=nmat,Nmat=Nmat,
                      ngroup=ngroup,nheight=nheight,
                      pi1=pi1,nloc=nloc)
        
        #store results after burn-in if they correspond to higher likelihood
        if (logl>max.llk & i>burnin){
          max.llk=logl
          MLE.pi1=pi1
          MLE.theta=theta
          MLE.z=z
          MLE.iter=i
        }
        
        #store results
        store.logl[i]=logl
        store.pi1[i,]=pi1
        store.theta[i,]=theta
    }
    
    #output MCMC results after discarding the burn-in phase
    #output MLE results
    seq1=burnin:ngibbs
    list(pi1=store.pi1[seq1,],
         theta=store.theta[seq1,],
         logl=store.logl,
         MLE.pi1=MLE.pi1,
         MLE.theta=MLE.theta,
         MLE.z=MLE.z,
         MLE.iter=MLE.iter,
         last.z=z)
}