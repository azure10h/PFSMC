#'Residual Resampling
#'
#'@description Also called remainder sampling, assuming that for particles with
#'large weights new particle can be assighed without drawing. In the
#'frist part of resampling step, samples with weights larger than 1/N
#'will be kept and reduced by a multiple of 1/N. In the second part,
#'we normalize the weights an do simple resampling.
#'@param w Sample weights
#'@param N Number of particles. Defalut is the length of w.
#'@return indx Return the chosen index.
#'@export
#'
resampleResidual<- function(w,N) {
  N=length(w)
  M=length(w)
  w=w/sum(w)
  indx=numeric(N)

  #integer parts
  Ns=floor(N*w)
  R=sum(Ns)
  #draw the determinstic part
  i=1
  j=0
  while(j<M){
    j = j + 1
    cnt=1
    while( cnt<= Ns[j]){
      indx[i]=j
      i = i+1
      cnt = cnt+1
    }
  }
  #Resample Multinomial
  N_rdn = N - R
  Ws=(N*w-Ns)/N_rdn
  Q=cumsum(Ws)
  while (i<=N) {
    sampl=runif(1)
    j=1
    while(Q[j]<sampl){
      j=j+1
    }
    indx[i]=j
    i=i+1
  }
  return(indx)
}

