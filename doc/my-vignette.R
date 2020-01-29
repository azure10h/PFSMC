## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
loss=function(theta,y) {
  return((theta-y)^2/2)
}

## -----------------------------------------------------------------------------
exploss = function(theta,y,mode=1,eta=1) {
  #mode = 1: first specify loss function then exp it
  #mode = 0: directly calculate the "likelihood"
  if(mode==1) {
    el=exp(-eta*loss(theta,y));
    } else {
      el=dnorm(y,theta,1);
  }
  return(el)
}

## -----------------------------------------------------------------------------
transition=function(theta,alpha) {
  n=length(theta) 
  u=runif(n)
  u=as.numeric(u<=(1-alpha))
  theta_new=theta*u+(a+(b-a)*runif(n))*(1-u)
  return(theta_new)
}

## -----------------------------------------------------------------------------
MH_move=function(theta,targetdist,prop_mean,prop_sig) {
  siz=length(theta)
  theta_new=rnorm(siz,prop_mean,prop_sig) #Generate samples from proposal distribution
  a1=targetdist(theta_new)/targetdist(theta) 
  a2=dnorm(theta,prop_mean,prop_sig)/dnorm(theta_new,prop_mean,prop_sig)
  a=a2*a1 
  #Accept the new candidate with probability min(1,a); otherwise just stay at the previous state.
  accept=rep(0,length(a)) 
  for (i in 1:length(a)) {
    accept[i]=min(1,a[i]) #alpha = min(1,a)
  }
  u=runif(siz) #The probability of accepting the new candidate
  u=as.numeric(u<accept)
  theta_new=theta_new*u+theta*(1-u)
  return(theta_new) #Return the Markov chain.
}

## -----------------------------------------------------------------------------
resampleMultinomial=function(weight) {
  M=length(weight) 
  Q=cumsum(weight)
  Q[M]=1
  indx=rep(0,M) #store sample index
  i=1
  while (i<=M) { 
    sampl=runif(1)
    j=1;
  #Resample according to sample weights
  while(Q[j]<sampl) {
    j=j+1
    }
  indx[i]=j
  i=i+1
  } 
  return(indx)
}

## ----echo=TRUE,eval=FALSE-----------------------------------------------------
#  devtools::install_github("azure10h/PFSMC")
#  library(PFSMC)

## -----------------------------------------------------------------------------
a=-10;b=10;N=1000;T=200 #Setting parameter space as (-10,10); particle numners N=1000; time=200
tha_sample=-a+(b-a)*runif(1000) #Initial particles theta_1^i for i=1:N
weight=matrix(0,T,N) #Store weights at each time
weight[1,]=rep(1,N)/N #Initial weights W_1^i=1/N

## -----------------------------------------------------------------------------
Z=rep(0,T) #Normalizing constants denoted in step 5
Z[1]=1 #Sum of weights is equal to 1 with equal weights
ESS=rep(0,T) #Store effctive sample sizes
ESS[1]=N #Initial effective sample size N
resample=numeric(T) #Store resample flags
theta_hat=numeric(T) #Store predicted parameters

## -----------------------------------------------------------------------------
datagenNoml = function(T,k,a,b,sig=1) { 
  if(k==0) {
    theta_true=(a+(b-a)*runif(1))*rep(1,T)
    Y=rnorm(T,theta_true,sig*rep(1,T))
  }
  #determine the chage point randomly
  change_point=floor(T/(k+1))*c(1:k)
  change_point=c(1,change_point,T)
  theta_list=a+(b-a)*runif(k+1)
  theta_true=rep(0,T)
  #Generate the true parameter
  for(j in 1:(k+1)) {
    theta_true[change_point[j]:change_point[j+1]]=theta_list[j]
  }
  #Generate data
  Y=rnorm(T,theta_true,sig*rep(1,T))
  return(list(Y,theta_true))
}

## ----fig.keep='none', message=FALSE, warning=FALSE,eval=FALSE-----------------
#  library(PFSMC)
#  set.seed(0421)
#  a=-10;b=10 #Parameter space
#  Data=datagenNoml(T,k,a,b)
#  Y=Data[[1]] #Data sequence
#  theta_true=Data[[2]] #True theta
#  Simulation<-PFSMC(Y=Y,eta=0.1,alpha=0.025,N=1000,c=0.5,T)
#  samples=Simulation$samples #Predictive distributions

## ----fig.keep='none',eval=FALSE-----------------------------------------------
#  ESS=Simulation[[1]] #Output: efective sample sizes.
#  Z=Simulation[[2]] #Output: normalizing constants.
#  theta_hat=Simulation[[3]] #Output: predictive parameters.
#  resampleflags=Simulation[[4]] #Output: resample flags.
#  plot(ESS,type = 'l',main = 'ESS',xlab = 'Time Index t')
#  abline(h=500,col='red',lwd=2)

## ----fig.keep="none",eval=FALSE-----------------------------------------------
#  plot(Y,type = "l",col="lightblue",axes=F,ylab="",xlab="")
#  par(new=T)
#  plot(data.frame(x=1:201,y=theta_true),type='l',col='red',ylim=c(a,b),
#       ylab='Data sequence/True mean/Predictive mean',xlab='time index')
#  title(main="Normal Data Sequence")
#  par(new=T)
#  plot(data.frame(x=1:201,y=theta_hat),type='l',col='blue',ylim=c(a,b),
#       ylab='',xlab='',axes=F)
#  ws=numeric(201)
#  for (i in 1:201){ws[i]=Simulation$samples[i,]%*%Simulation$weight[i,]}
#  par(new=T)
#  plot(data.frame(x=1:201,y=ws),type='l',col='green',
#       ylim=c(-10,10),ylab='',xlab='',axes=F)
#  legend("bottomright",legend = c("Data","True mean","Predictive mean","Weighted mean"),lty=c(1,1,1),
#         col=c("lightblue","blue","red","green"),cex = 0.5)

