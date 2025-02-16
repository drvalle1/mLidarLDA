#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
// Sample from categorical distribution

//Input:
// value this is the value of a uniform random variable 
// prob vector of the probabilities associated with each category

//Output:
// res selected category
int whichLessDVPresenceFast(double value, 
                            NumericVector prob) {
  int res=-1;
  double probcum = 0;
  
  for (int i = 0; i < prob.length(); i++) {
    probcum = probcum + prob[i];
    if (value < probcum) {
      res = i;
      break;
    }
  }
  return res;
}

// [[Rcpp::export]]
// This function samples the z vector containing the FST allocated to each pixel

//Inputs:
// ltheta log of the J x K theta matrix containing the proportion of each FST k in each transect j
// nmat P x H matrix containing the number of returns for each pixel p in each height class h
// Nminusn P x H matrix calculated as Nmat-nmat
// lpi1 log of the K x H pi matrix containing the return probability of each FST k for each height class h
// l1minuspi1 K x H matrix calculated as log of 1-pi1
// ngroup maximum number of FSTs
// nloc number of pixels
// nheight number of height classes
// randu vector of uniform random variables
// TransID transect identifier determining which pixel belong to which transect

//Output:
//zvec vector containing the FST allocated to each pixel

IntegerVector samplez(NumericMatrix ltheta,
                      IntegerMatrix nmat,
                      IntegerMatrix Nminusn,
                      NumericMatrix lpi1,
                      NumericMatrix l1minuspi1,
                      int ngroup,
                      int nloc,
                      int nheight,
                      NumericVector randu,
                      IntegerVector TransID) {

  IntegerVector zvec(nloc);
  NumericVector prob(ngroup);
  int znew;
  double soma=0;
   
  for(int i=0; i<nloc;i++){
    for (int k=0; k<ngroup; k++){
      soma=0;
      for (int j=0; j<nheight; j++){
        soma=soma+nmat(i,j)*lpi1(k,j)+
                  Nminusn(i,j)*l1minuspi1(k,j);
      }
      prob[k]=soma+ltheta(TransID[i],k);
    }
    prob=prob-max(prob);
    prob=exp(prob);
    prob=prob/sum(prob);

    //multinomial draw
    znew=whichLessDVPresenceFast(randu[i],prob);
    zvec[i]=znew+1;
  }
  return zvec;
}

// [[Rcpp::export]]
// This function tabulates the number of pixels in each transect assigned to each FST

//Inputs:
// z vector containing the FST allocated to each pixel
// ngroup maximum number of FSTs
// nloc number of pixels
// ntransect number of transects
// TransID transect identifier determining which pixel belong to which transect

//Output:
//ntk table containing the number of pixels in each transect assigned to each FST

IntegerMatrix calc_ntk(IntegerVector z,
                      int ngroup,
                      int nloc,
                      int ntransect,
                      IntegerVector TransID) {
  
  IntegerMatrix ntk(ntransect,ngroup);

  for(int i=0; i<nloc;i++){
    ntk(TransID[i],z[i])=ntk(TransID[i],z[i])+1;
  }
  return ntk;
}