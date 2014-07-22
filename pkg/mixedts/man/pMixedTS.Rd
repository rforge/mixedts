\name{pMixedTS}
\alias{pMixedTS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ MixedTS Cumulative distribution function.
%%  ~~function to do ... ~~
}
\description{Cumulative Distribution function of mixed Tempered Stable distribution
}
\usage{
pMixedTS(q, mu0, mu, sig, a, alpha, lambda_p, lambda_m,
         lower = -Inf, setSup = NULL, setInf = NULL, MixingDens = "Gamma")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{q}{ vector of quantiles.
%%     ~~Describe \code{q} here~~
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
  \item{alpha}{  Density parameter ranging from 0 to 2.
%%     ~~Describe \code{alpha} here~~
}
  \item{lambda_p}{ Tempering parameter for the right tail, its value is graeter than 0.
%%     ~~Describe \code{lambda_p} here~~
}
  \item{lambda_m}{ Tempering parameter for the right tail, its value is graeter than 0.
%%     ~~Describe \code{lambda_m} here~~
}
  \item{lower}{ Internal parameter.
%%     ~~Describe \code{lower} here~~
}
  \item{setSup}{ Internal parameter.
%%     ~~Describe \code{setSup} here~~
}
  \item{setInf}{ Internal parameter.
%%     ~~Describe \code{setInf} here~~
}
  \item{MixingDens}{ Internal parameter.
%%     ~~Describe \code{MixingDens} here~~
}
}
\details{ This function returns the cdf of mixed tempered stable random variable.
%%  ~~ If necessary, more details than the description above ~~
}
%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%}
%\references{
%% ~put references to the literature/web site here ~
%}
%\author{
%%  ~~who you are~~
%}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
x<-seq(-3,3,length=400) 
prob <-pMixedTS(x,mu0=-0.077,mu=0.027,sig=1.015,a=1.269,alpha=1.372,
                lambda_p=1.025,lambda_m=0.996,lower=-Inf,
                            setSup=NULL,setInf=NULL,MixingDens="Gamma")
plot(x,prob,type="l")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
