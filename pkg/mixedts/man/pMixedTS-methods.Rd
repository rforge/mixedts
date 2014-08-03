\name{pMixedTS-methods}
\docType{methods}
\alias{pMixedTS}
\alias{pMixedTS-methods}
\alias{pMixedTS,param.MixedTS-method}
\title{Probability of Mixed Tempered Stable distribution}
\description{
This Method returns the cdf of a Mixed Tempered Stable
}
\section{Methods}{
\describe{

\item{\code{signature(object = "param.MixedTS",x = numeric(), setSup=NULL,setInf=NULL,N=2^10)}}{
This method returns an object of class \code{MixedTS} where the slot \code{prob} contains the value of the probability evaluated on the \code{x}. \code{setSup} and \code{setInf} are used to choose \code{+ infinity} and \code{- infinty}. \code{N} is the number of point used for discretization in \code{fft} algorithm.
}
}}


\examples{
# First Example

# Density of MixedTS with Gamma

ParamEx1<-setMixedTS.param(mu0=0, mu=0, sigma=0.4, a=1.5,
                           alpha=0.8, lambda_p=4, lambda_m=1, 
                           Mixing="Gamma")
# support

x<-seq(-3,1,length=100)

prob1<-pMixedTS(x=x,object=ParamEx1,setSup=10,setInf=-10,N=2^7)

plot(prob1)

# Prob of MixedTS with IG

Mix<-"User"

parMix<-list(lamb=1,mu1=1)

logmgf<-("lamb/mu1*(1-sqrt(1-2*mu1^2/lamb*u))")

ParamEx2<-setMixedTS.param(mu0=0, mu=0, sigma=0.4, a=logmgf,
                           alpha=0.8, lambda_p=4, lambda_m=1,
                           Mixing=Mix,paramMixing=parMix)

x<-seq(-3,1,length=100)

prob2<-pMixedTS(x=x,object=ParamEx2,setSup=10,setInf=-10,N=2^7)
plot(prob2)


}


