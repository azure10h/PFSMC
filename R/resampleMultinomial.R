#'A multinomial resampling function
#'@description To construct an unweighted sample from a weighted one.
#' Resample will be constructed according to the samples' weight.
#' Samples with higher weights will be selected more often than
#' samples with lower weights.
#'@param w weights vector.
#'@return Function return an index called ancestry vector, which is
#' drawn from a categorical distribution.
#' Replace samples with the ancestor.
#'@export



resampleMultinomial=function(w){

  M=length(w)
  Q=cumsum(w)
  Q[M]=1
  indx=rep(0,M)

  i=1
  while (i<=M) { sampl=runif(1)
  j=1;
  while(Q[j]<sampl){j=j+1}
  indx[i]=j
  i=i+1
  }
  return(indx)
}
