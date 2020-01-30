#'Residual Resampling
#'@export
resampleResidual<- function(w,N)
{
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

# w=runif(100)
# N=100
# resampleResidual(w,N)
