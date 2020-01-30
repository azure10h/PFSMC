#'Sequential Monte Carlo Partical Filter
#'@description A kinetic prediction implemented
#'with Sequential Monte Carlo algorithm
#'
#'
#'@param Y a 1*T vector of sequential data sequence.
#'@param eta learning parameter. Determines the rate of weight updating process.
#'@param alpha mixing parameter. Determines the speed of model convergence and the rate of trakcing changepoints.
#'@param N number of particles to predict underlying densities.
#'@param c effective sample size thershold.
#'@param T number of data sequence. Time index.
#'@param loss loss function for the underlying density. Used to update sample weights.
#'@param resample resampling function. 
#'
#'@return \code{PFSMC} returns a list of effective sample size, normalized constants, predicted parameters theta and resample flags at each time.
#'@export
#'
#'@examples
#'
#'
#'#Generate true parameters and data sequence with 5 change points
#'k=5; T=200; a=-10; b=10
#'library(PFSMC)
#'Data=datagenNoml(200,5,-10,10)
#'Y=Data[[1]]
#'theta_true=Data[[2]]
#'
#'#Detecting changepoints using `PFSMC` funciton.
#'
#'Simulation<-PFSMC(Y=Y,eta=10*T^(-1/3),alpha=k/(T-1),N=1000,
#'c=0.5,T=200,loss=lossGaussian,resample=resampleMultinomial)
#'ESS=Simulation[1]
#'theta_hat=Simulation$theta_hat
#'
#'#Result visulization
#'plot(theta_true,type="l",ylim=c(-12,12))
#'lines(theta_hat,col="red")



PFSMC=function(Y,eta,alpha,N,c,T,loss,resample=resampleMultinomial)
{

  loss=match.fun(loss)
  resampling=match.fun(resample)

  mode=1
  samples=matrix(0,T,N)
  tha_sample=a+(b-a)*runif(N)    #Initial particles theta_1^i for i=1:N
  samples[1,]=tha_sample         #store particle samples
  weight=matrix(0,T,N)           #Initial weights W_1^i=1/N
  weight[1,]=rep(1,N)/N

  Z=rep(0,T) #Normalizing constants denoted in step 5
  Z[1]=1 #Sum of weights is equal to 1 with equal weights
  ESS=rep(0,T) #Store effctive sample sizes
  ESS[1]=N #Initial effective sample size N
  resample=numeric(T) #Store resample flags
  theta_hat=numeric(T) #Store predicted parameters
  theta_hat[1]=mean(tha_sample) #Initial predictive parameter


#filteringdst_part1 function
#This function is important in MH calculating step
filteringdst_part1=function(t,alpha,Z)
  {
  Coef=matrix(0,t,1)
  Coef[1]=(1-alpha)^(t-1)
  if(t>=2)
    {
    for(k in 2:t)
    {
      Coef[k]=alpha*(1-alpha)^(t-k)*Z[k]
    }
  }
return(Coef)
  }

# filteringdst_part2 function
filteringdst_part2=function(theta,t,Y,mode=1,eta=1)
{

  if(is.null(nrow(theta))) {theta=t(t(theta))}
  X=matrix(0,length(theta),t)

  if(mode==0){
    X[,t]=exploss(theta,Y[t],mode,eta)
    k=t-1
    while(k>=1){
      X[,k]=X[,k+1]*exploss(theta,Y[k],mode,eta)
      k=k-1
    }
  }
  else {
      X[,t]=loss(theta,Y[t])
      k=t-1
      while(k>=1){
        X[,k]=X[,k+1]+loss(theta,Y[k])
        k=k-1
      }
      X=exp(-eta*X)
    }
  return(X)
}


#MH_move function
#This function realize the moving step using the Metropolis-Hastings algorithm.
#We view the sampling process as a hidden state which follows a MCMC kernel.
MH_move=function(theta,targetdist,prop_mean,prop_sig) {
  siz=length(theta)
  theta_new=rnorm(siz,prop_mean,prop_sig) #Generate samples from a proposal distribution
  a1=targetdist(theta_new)/targetdist(theta) #a1 is the target distribution we want to sample from
  a2=dnorm(theta,prop_mean,prop_sig)/dnorm(theta_new,prop_mean,prop_sig) #a2 is the values we sample from the proposal distribution
  a=a2*a1
  #Accept the new candidate with probability min(1,a); otherwise just stay at the previous state.
  accept=rep(0,length(a))
  for (i in 1:length(a)) {
    accept[i]=min(1,a[i]) #alpha = min(1,a)
  }
  u=runif(siz) #The probability of accepting the new candidate
  u=as.numeric(u<accept)
  theta_new=theta_new*u+theta*(1-u) #Move the samples
  return(theta_new) #Return the Markov chain.
}

#transition function
#We use the transition function to achieve the mixing step.
transition=function(theta,alpha,a,b) {
  n=length(theta)
  u=runif(n)
  #With probability alpha, we keep the new particle; with probability alpha,
  #replace it with a random draw from parameter space.
  u=as.numeric(u<=(1-alpha))
  theta_new=theta*u+(a+(b-a)*runif(n))*(1-u)
  return(theta_new)
}



#exploss function
#Based on the loss funciton we input, calculate the exponential loss of samples.
exploss = function(theta,y,mode=1,eta=1)
{
  #\exp(-\eta loss(theta,y))
  #mode = 1: first specify loss function then exp it
  #mode = 0: directly calculate the "likelihood"

  if(mode==1) {
    el=exp(-eta*loss(theta,y));
  }
  else {el=dnorm(y,theta,1);}
  return(el)
}


isrejuvenate = 1
ismixing = 1

for (t in 1:(T-1)) {
  
  #update weight: particle with larger score gains larger weight
  weight[t+1,]=weight[t,]*(exploss(tha_sample,Y[t],mode,eta))
  Z[t+1]=sum(weight[t+1,])
  weight[t+1,]=weight[t+1,]/Z[t+1] #normalize the weights
  Z[t+1]=Z[t]*Z[t+1] #a useful constatnt, product sum of weights that will be used in calculating the integrals

  #calculate ESS
  ESS[t+1]=1/sum(weight[t+1,]^2)
  resample_flag=0

  if(ESS[t+1]<c*N) { 
    
    #resample when ESS fall below a certain threshold
    resample_flag=1  
    ind=resampling(weight[t+1,]) #Using a chosen resampling method to resample particles.
    tha_sample=tha_sample[ind]  #Decide which samples need to be resampled.
    weight[t+1,]=rep(1,N)/N   #Obtain equal weights.

    #rejuvenate/Move using MH kernal
    if(isrejuvenate) {
    
      prop_mean=mean(tha_sample) 
      prop_sig=sd(tha_sample)
      Coef=filteringdst_part1(t,alpha,Z)
      Xcal=function(theta) {
        filteringdst_part2(theta,t,Y,mode,eta)
        }
      targetdist=function(theta) {
          return(t(Xcal(theta)%*%Coef))
          }
      tha_sample=MH_move(tha_sample,targetdist,prop_mean,prop_sig)
    }
  }

  #mixing step: move samples according to a mixing transition kernel
  if(ismixing) {
    tha_sample=transition(tha_sample,alpha,a,b)
  }

  if(resample_flag) {
    theta_hat[t+1]=mean(tha_sample) #we use the mean of samples as predictive value
    samples[t+1,]=tha_sample 
  }
  else {
    theta_hat[t+1]=mean(tha_sample[resampling(weight[t+1,])])
    samples[t+1,]=tha_sample[resampling(weight[t+1,])]
    }
  
  resample[t+1]=resample_flag #record resample flags
}
return(list(ESS=ESS,Z=Z,theta_hat=theta_hat,resample_flag=resample,samples=samples,weight=weight))
}

