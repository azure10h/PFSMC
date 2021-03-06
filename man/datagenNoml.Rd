% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/datagennoml.R
\name{datagenNoml}
\alias{datagenNoml}
\title{Data Generating Function}
\usage{
datagenNoml(T, k, a, b, sig = 1)
}
\arguments{
\item{T}{Index of time series data.}

\item{k}{Number of abrupt changes}

\item{a}{Parameter space, upper bound.}

\item{b}{Parameter space, lower bound.}

\item{sig}{Variance of normal distribution with 1 as default.}
}
\value{
Return `theta_true` as ture underlying paramaeters and
`Y` as data sequence.
}
\description{
A function creating parameters with abrupt
changes and sequential data generated by those parameters
in Gaussian form.
}
\examples{

data=datagenNoml(T=1000,k=6,a=-10,b=10,sig=1)
y=data[[1]]
theta_t=data[[2]]
plot(data.frame(x=1:T,y=y),type='l',col='lightblue',
ylim=c(-10,10),ylab='Generate Data/True mean',xlab='time index')
par(new=TRUE)
plot(data.frame(x=1:T,y=theta_t),type='l',
col='red',ylim=c(-10,10),ylab='',xlab='',axes=F)
title(main="Normal Data Sequence")

}
