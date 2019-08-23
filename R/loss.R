#'Loss Function
#'
#'@param theta Underlying parametes theta. 
#'@param y Data sequence generated from underlying distributions.
#'@description A loss function calculating the squared score of parameters and data sequence.
#'@return Returns a squared loss. 
#'@export
#'@examples
#'
#'#First generate data sequence with T=100, parameter space=[-10,10]
#'library(PFSMC)
#'data=datagenNoml(T=100,k=6,a=-10,b=10,sig=1)
#'y=data[[1]]; theta_t=data[[2]]
#'loss1=numeric(100)
#'for (i in 1:100)
#'{ loss1[i]=loss(theta_t[i],y[i]) }
#'loss1
#'
#'
#'
loss=function(theta,y)
{
  return((theta-y)^2)
}
