\name{mle.MixedTS}
\alias{mle.MixedTS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{MLE estimation
%%  ~~function to do ... ~~
}
\description{The code performs the maximum likelihood estimation.
%%  ~~ A concise (1-5 ines) description of what the function does. ~~
}
\usage{
mle.MixedTS(data, param0, method = "BFGS", fixed.param, 
            lower.param, upper.param, MixingDens = "Gamma",N=2^7)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ vector contains the dataset.
%%     ~~Describe \code{data} here~~
}
  \item{param0}{Starting point for optimization routine
%%     ~~Describe \code{param0} here~~
}
  \item{method}{ Optimization Method
%%     ~~Describe \code{method} here~~
}
  \item{fixed.param}{ parameters that the user chooses to fix during the optimization routine.
%%     ~~Describe \code{fixed.param} here~~
}
  \item{lower.param}{ Lower bounds of parameters.
%%     ~~Describe \code{lower.param} here~~
}
  \item{upper.param}{ Upper bounds for parameters.
%%     ~~Describe \code{upper.param} here~~
}
  \item{MixingDens}{ Internal parameter.
%%     ~~Describe \code{MixingDens} here~~
}
\item{N}{ Internal parameter.
%%     ~~Describe \code{MixingDens} here~~
}
}
\details{ The function returns the estimated parameters using the maximum likelihood estimation procedure
%%  ~~ If necessary, more details than the description above ~~
}
\value{ Value of Estimated Parameters, Variance-Covariance Matrix and log-likelihood.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
# generation of random numbers.
set.seed(123)
NumberMixedTS <- rMixedTS(n=1000,
                          mu0=-0.077,mu=0.027,sig=1.015,a=1.269,
                          alpha=1.372,lambda_p=1.025,lambda_m=0.996,
                          setSup=NULL,setInf=NULL,MixingDens="Gamma")
# Starting point
param0<-c(0,0,sig=0.9,a=1,
              alpha=1.05,lambda_p=1.1,lambda_m=0.9)

names(param0) <- c("mu0","mu","sig","a","alpha","lambda_p","lambda_m")

#Maximum Likelihood Estimation

timeMixedTS<-system.time(
result<-mle.MixedTS(NumberMixedTS,param0=param0, 
                     fixed.param=NULL,
                    lower.param=NULL,
                    upper.param=NULL,
                    MixingDens="Gamma",N=2^7)
)

# Results
result$estLevpar

sqrt(diag(result$covLev))

result$logLik

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
