\name{rMixedTS}
\alias{rMixedTS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Random Number
%%  ~~function to do ... ~~
}
\description{This function returns a set of random numbers generated from a MixedTS distribution.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
rMixedTS(n, mu0, mu, sig, a, alpha, lambda_p, lambda_m, 
         setSup = NULL, setInf = NULL, MixingDens = "Gamma")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{ length of the vector that contains the generated random number.
%%     ~~Describe \code{n} here~~
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
 
}
\value{Vector containing random numbers of a MixedTS distribution.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}

\examples{
set.seed(123)
NumberMixedTS <- rMixedTS(n=1000,
                          mu0=-0.077,mu=0.027,sig=1.015,a=1.269,
                          alpha=1.372,lambda_p=1.025,lambda_m=0.996,
                          setSup=NULL,setInf=NULL,MixingDens="Gamma")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
