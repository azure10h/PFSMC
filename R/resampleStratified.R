#'Stratified Resampling
#'@description A resampling function. Stratified resampling is used to
#'reduce variance when Monte Carlo methods are implemented to estimate
#'population statistics. We pre-partition the (0,1] into n disjoint sets,
#'then draw independetly from each sub-intervals.
#'
#'@param w Sample weights.
#'@param N Number of particles. Defalut is the length of w.
#'@return Return the index that is chosen.
#'
#'
#'@export

resampleStratified=function(w,N){

  N=length(w)
  M=length(w)
  w=w/sum(w)
  Q=cumsum(w)
  indx=numeric(N)
  T=seq(0,1-1/N,length.out = N)+runif(N)/N

  i = 1
  j = 1
  while(i<=N && j<=M){
    while(Q[j]<T[i]){
      j=j+1
    }
    indx[i]=j
    i = i+1
  }
  return(indx)
}
