#'Systematic Resampling
#'
#'@description Similat to the stratified resampling. Difference is,
#'through whole resampling steps, random number is drawed only once.
#'Refer to \code{\link{resampleStratified}}
#'
#'@param w Sample weights
#'@param N Number of particles. Defalut is the length of w.
#'
#'@return index Return the chosen index of samples.
#'@export

resampleSystematic=function(w,N) {
  N=length(w)
  M=length(w)
  w=w/sum(w)
  Q=cumsum(w)
  indx=numeric(N)
  T=seq(0,1-1/N,length.out = N)+runif(1)/N

  i=1
  j=1

  while(i<=N && j<=M) {
    while (Q[j]<T[i]) {
      j=j+1
    }
    indx[i]=j
    i=i+1
  }
return(indx)
  }


