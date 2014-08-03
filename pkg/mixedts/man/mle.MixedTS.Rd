\name{mle.MixedTS}
\alias{Mixed Tempered Stable distribution}
\alias{Normal Variance Mean Mixture}
\alias{mle.MixedTS}
\alias{mle}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Maximum Likelihood Estimation for MixedTS distribution
}
\description{
Estimate MixedTS parameters using the Maximum Likelihood Estimation procedure.
}
\usage{
mle.MixedTS(object, start = list(), Data = NULL, 
          method = "L-BFGS-B", fixed.param = NULL, 
          lower.param = NULL, upper.param = NULL, 
          setSup = NULL, setInf = NULL, N = 2^10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{ 
  \item{object}{ an object of class \code{param.MixedTS} that contains informations about the model. 
%%     ~~Describe \code{object} here~~
}
  \item{start}{ a list of parameter for the mle.
%%     ~~Describe \code{start} here~~
}
  \item{Data}{ a numeric object containing the dataset.
%%     ~~Describe \code{Data} here~~
}
  \item{method}{ methods for optimization routine. See \code{optim} for more details. 
%%     ~~Describe \code{method} here~~
}
  \item{fixed.param}{ a list of the model parameter that must be fix during optimization routine. Choosing \code{alpha=2} the function returns the estimate parameters for the Normal Variance Mean Mixture distribution.
%%     ~~Describe \code{fixed.param} here~~
}
  \item{lower.param}{ a list containing the lower bound for the parameters.
%%     ~~Describe \code{lower.param} here~~
}
  \item{upper.param}{ a list containing the upper bound for the parameters.
%%     ~~Describe \code{upper.param} here~~
}
  \item{setSup}{ Internal parameter. see documentation for \code{dMixedTS} for more details.
%%     ~~Describe \code{setSup} here~~
}
  \item{setInf}{ Internal parameter. see documentation for \code{dMixedTS} for more details.
%%     ~~Describe \code{setInf} here~~
}
  \item{N}{ Internal parameter. see documentation for \code{dMixedTS} for more details.
%%     ~~Describe \code{N} here~~
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
The function returns an object of class \code{MixedTS.qmle}.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\examples{
# First Example:
# We define the Mixed Tempered Stable using the function setMixedTS.param


ParamEx1<-setMixedTS.param(mu0=0, mu=0, sigma=0.4, a=1.5,
                           alpha=0.8, lambda_p=4, lambda_m=1, Mixing="Gamma")

# We generate a sample using the rMixedTS method
set.seed(100)
Rand1 <- rMixedTS(x=5000,object=ParamEx1, setSup=10,setInf=-10,N=2^9)

# Estimate procedure
\dontrun{
est1<-mle.MixedTS(object=Rand1 , setSup=10,setInf=-10,N=2^9)
# Show results

summary(est1)
}



}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ M }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
