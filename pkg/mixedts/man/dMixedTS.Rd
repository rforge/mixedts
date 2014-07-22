\name{dMixedTS}
\alias{dMixedTS}
\alias{density of Mixed Tempered Stable}
\alias{Mixed TS density.}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Mixed Tempered Stable Density function
%%  ~~function to do ... ~~
}
\description{Density function of mixed Tempered Stable distribution.
}
\usage{
dMixedTS(x, mu0, mu, sig, a, alpha, lambda_p, lambda_m, 
         setSup = NULL, setInf = NULL, MixingDens = "Gamma",N=2^10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ vector of quantiles.
%%     ~~Describe \code{x} here~~
}
  \item{mu0}{ Density parameter, It ranges from -infinity to +infinity.
%%     ~~Describe \code{mu0} here~~
}
  \item{mu}{ Density parameter ranging from -infinity to +infinity. 
%%     ~~Describe \code{mu} here~~
}
  \item{sig}{ Density parameter ranging from 0 to +infinity.
%%     ~~Describe \code{sig} here~~
}
  \item{a}{ shape parameter of the mixing density (greater than zero).
%%     ~~Describe \code{a} here~~
}
  \item{alpha}{ Density parameter ranging from 0 to 2.
%%     ~~Describe \code{alpha} here~~
}
  \item{lambda_p}{ Tempering parameter for the right tail, its value is graeter than 0.
}
  \item{lambda_m}{ Tempering parameter for the right tail, its value is graeter than 0.
}
  \item{setSup}{ Internal parameter.
}
  \item{setInf}{ Internal parameter.
}
  \item{MixingDens}{Internal parameter. 
}
  \item{N}{Internal parameter. 
}
}
\details{This function returns the density of mixed tempered stable random variable.
%%  ~~ If necessary, more details than the description above ~~
}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
# code for testing the density of Mixed TS
x<-seq(-3,3,length=400) 
y<-dMixedTS(x,mu0=-0.077,mu=0.027,sig=1.015,a=1.269,alpha=1.372,lambda_p=1.025,lambda_m=0.996)

plot(x,y,type="l")

}
