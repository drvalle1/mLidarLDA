
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Author

Denis Valle

# Modified LidarLDA

<!-- badges: start -->
<!-- badges: end -->

The goal of the modified LidarLDA is to identify Forest Structural Types
(FSTs), estimating two sets of parameters: one set characterizes each
transect in relation to the relative abundance of FSTs (`theta` matrix)
while the other set characterizes each FST in relation to its vertical
structure based on the return probability for each height stratum (`pi1`
matrix).

In this tutorial, we describe how this model can be fit and the
corresponding results displayed and interpreted using simulated data.

## Simulated data

The simulated dataset contains information for 4,000 pixels evenly
distributed in 100 LiDAR transects. Furthermore, there are 10 height
bins, labeled h1, h2, …, h10. These data were simulated with only 2 FSTs
and the code used to generate these data is in the R script ‘simulated
data.R’.

The data in `simdat nmat1.csv` consist of the number of returned light
pulses in each pixel (rows; p1,…,pP) and each height class (columns;
h1,…,hH).

The data in `simdat NMAT.csv` consist of the corresponding number of
incoming light pulses in each pixel (rows; p1,…,pP) and each height
class (columns; h1,…,hH). As a result, the values in `simdat nmat1.csv`
are always smaller or equal to those in `simdat NMAT.csv`.

``` r
rm(list=ls()) #clean the environment

#get simulated data
setwd('/Users/denisvalle/Documents/EBA stuff/simulated example')
nmat=read.csv('simdat nmat1.csv')
head(nmat)
#>    X h1 h2 h3 h4 h5 h6 h7 h8 h9 h10
#> 1 p1 55 88 77 68 48 47 55 44 45  34
#> 2 p2 64 69 88 74 51 58 47 60 34  35
#> 3 p3 77 88 78 85 68 60 48 54 31  37
#> 4 p4 54 80 79 67 64 61 56 47 37  24
#> 5 p5 48 90 69 67 60 51 50 46 29  35
#> 6 p6 54 87 68 65 60 42 55 36 34  20

Nmat=read.csv('simdat NMAT.csv')
head(Nmat)
#>    X  h1  h2  h3  h4  h5  h6  h7  h8  h9 h10
#> 1 p1  93 102 109 100  87  85 103  98 107 109
#> 2 p2 113  87 117  94  81  96 100 111  91 105
#> 3 p3 112 103 107 116 109 106  97 103  97 109
#> 4 p4 104  95 112 105 104 107 110  92  78  79
#> 5 p5  84 112  90  99 104  94 103 108  85  93
#> 6 p6 104 106 101 106 105  87 104 102  89  81

npix=nrow(nmat); npix #number of pixels
#> [1] 40000
nheight=ncol(nmat)-1; nheight#number of height classes
#> [1] 10

#eliminate row names contained in first column and 
#convert to data matrix
nmat=data.matrix(nmat[,-1])
Nmat=data.matrix(Nmat[,-1])
```

We also import an auxiliary file called `simdat auxiliary.csv`. This
file contains information regarding which pixels (`pixel.id`) belong to
which transects (`trans.id`) as well as the spatial coordinates of each
transect (`x` and `y`).

``` r
#get auxiliary file
setwd('/Users/denisvalle/Documents/EBA stuff/simulated example')
aux=read.csv('simdat auxiliary.csv')
head(aux)
#>   trans.id x y pixel.id
#> 1        1 1 1        1
#> 2        1 1 1        2
#> 3        1 1 1        3
#> 4        1 1 1        4
#> 5        1 1 1        5
#> 6        1 1 1        6
```

## Installation

The code described in this tutorial requires the following R packages:

- `readr`
- `ggplot2`
- `Rcpp`

``` r
library('readr')
library('ggplot2')
library('Rcpp')
```

You can obtain the required functions for our algorithm by downloading
the files in the public repository
[mLidarLDA](https://github.com/drvalle1/mLidarLDA). Once these files
have been downloaded, they can be read using the commands below within
R:

``` r
setwd('/Volumes/Users/drvalle/GIT_models/mLidarLDA')
sourceCpp('code.cpp')
source('mLidarLDA_main.R')
source('mLidarLDA_functions.R')
```

## Fitting the model

We fit these simulated data using the code below. We rely on the
following settings:

- `ngroup=5`: we assume a maximum of 5 FSTs knowing that the model can
  identify fewer FSTs if the data provide support for this;
- `ngibbs=1000`: we set to 1000 the number of iterations of the gibbs
  sampler;
- `burnin=ngibbs*0.5`: the number of iterations to discard as burn-in
  was set to 500;
- `gamma1=0.1`: this parameter controls how strongly the model will seek
  to identify fewer groups. Values closer to 0 tend to correspond to
  fewer groups.

We fit the model using the function `mLidarLDA`. We export the following
results from the model:

- `logl`: this vector shows the log-likelihood at each iteration of the
  Gibbs sampler. This is useful to assess algorithm convergence;
- `MLE.pi1`: this is the `pi1` matrix associated with the highest
  likelihood after burn-in. The `pi1` matrix contains the return
  probability for each FST in each height class;
- `MLE.theta`: this is the `theta` matrix associated with the highest
  likelihood after burn-in. The `theta` matrix contains the proportion
  of each FST in each pixel.

``` r
#settings
gamma1=0.1
ngroup=5
ngibbs=1000
burnin=ngibbs*0.5

#fit model
mod=mLidarLDA(nmat=nmat,Nmat=Nmat,ngroup=ngroup,
              ntransect=ntransect,trans.id=trans.id,
              ngibbs=ngibbs,burnin=burnin,gamma1=gamma1)

#export results
setwd('/Users/denisvalle/Documents/EBA stuff/simulated example')
write.csv(mod$logl,'res logl.csv',row.names=F)
write.csv(mod$MLE.pi1,'res MLE pi1.csv',row.names=F)
write.csv(mod$MLE.theta,'res MLE theta.csv',row.names=F)
```

Next, we import these results to examine them more closely. For example,
we assess convergence by examining the trace-plot of the log-likelihood.
These plots suggest that the algorithm has converged after a few
iterations.

``` r
#get model results
setwd('/Users/denisvalle/Documents/EBA stuff/simulated example')
pi1.estim=data.matrix(read.csv('res MLE pi1.csv'))
theta.estim=data.matrix(read.csv('res MLE theta.csv'))
logl=read.csv('res logl.csv')

#get true parameter values
pi1.true=read.csv('simdat pi1.csv')
theta.true=read.csv('simdat theta.csv')

#assess convergence
par(mfrow=c(1,2))
plot(1:nrow(logl),logl[,1],type='l',xlab='Iterations',ylab='log-likelihood')
plot(200:nrow(logl),logl[200:nrow(logl),1],type='l',xlab='Iterations',ylab='log-likelihood')
```

<img src="man/figures/README-unnamed-chunk-7-1.png" width="100%" />

Furthermore, according to the estimated `theta` matrix, our model has
identified 2 (out of a maximum of 5) FSTs. These 2 FSTs, on average,
represent 99.99% of the proportion of FSTs.

``` r
colnames(theta.estim)=paste0('FST ',1:5)
boxplot(theta.estim,ylab='Proportion of FST')
```

<img src="man/figures/README-unnamed-chunk-8-1.png" width="100%" />

``` r
mean(apply(theta.estim[,1:2],1,sum))
#> [1] 0.9998664
```

Because this is based on simulated data, there is a nice pattern
regarding how the relative abundance of each FST changes in space. More
specifically, we simulated the data assuming that there is a river that
crosses the middle of the landscape horizontally and that FST 2 is much
more common along the river while FST 1 is more common further away from
the river. We separately display the map depicting the proportion of
each FST.

``` r
#get spatial coordinates for each transect
aux1=unique(aux[,c('x','y','trans.id')])
colnames(theta.estim)=paste0('FST',1:ncol(theta.estim))
theta.estim1=cbind(aux1,theta.estim)

#display estimated theta matrix for FST1
ggplot()+
  geom_tile(data=theta.estim1,alpha=0.5,aes(x=x,y=y,fill=FST1))+
  scale_fill_gradient2(low = "cyan", mid = "red",high='purple',
                       limits=c(0,1),midpoint=0.5) +
  ggtitle('FST 1')+
  theme(legend.title=element_blank())
```

<img src="man/figures/README-unnamed-chunk-9-1.png" width="100%" />

``` r

#display estimated theta matrix for FST2
ggplot()+
  geom_tile(data=theta.estim1,alpha=0.5,aes(x=x,y=y,fill=FST2))+
  scale_fill_gradient2(low = "cyan", mid = "red",high='purple',
                       limits=c(0,1),midpoint=0.5) +
  ggtitle('FST 2')+
  theme(legend.title=element_blank())
```

<img src="man/figures/README-unnamed-chunk-9-2.png" width="100%" />

Finally, we show below the estimated return probabilities for each FST
in each height class. Notice that, as expected for the vegetation close
to the river, FST 2 tends to have lower probabilities across most of the
height classes, except for the first height class. On the other hand,
FST 1 generally has higher return probabilities except for the first
class.

``` r
rownames(pi1.estim)=paste0('FST ',1:5)

#display estimated vertical profiles of each FST
plot(1:nheight,pi1.estim['FST 1',],type='h',col='red',ylim=c(0,1),
     xlab='height class (m)',ylab='Return probabilities',lwd=2)
for (i in 1:nheight){
  lines(rep(i,2)+0.1,c(0,pi1.estim['FST 2',i]),col='blue',lwd=2)
}

#add legend
legend(x=8,y=1,lty=1,lwd=2,legend=paste0('FST ',1:2),
       col=c('red','blue'))
```

<img src="man/figures/README-unnamed-chunk-10-1.png" width="100%" />

Finally, a comparison between the estimated and true `theta` and `pi1`
matrices is provided below. As expected, all of the parameters were well
estimated, with all points near the one-to-one diagonal red line.

``` r
rango=c(0,1)

#compare estimated vs true theta matrices
plot(data.matrix(theta.true[,-1]),
     data.matrix(theta.estim[,1:2]),xlim=rango,ylim=rango,
     main='theta matrix',xlab='true',ylab='estimated')
lines(rango,rango,col='red')
```

<img src="man/figures/README-unnamed-chunk-11-1.png" width="100%" />

``` r

#compare estimated vs true pi1 matrices
plot(data.matrix(pi1.true[,-1]),
     data.matrix(pi1.estim[1:2,]),xlim=rango,ylim=rango,
     main='pi1 matrix',xlab='true',ylab='estimated')
lines(rango,rango,col='red')
```

<img src="man/figures/README-unnamed-chunk-11-2.png" width="100%" />
