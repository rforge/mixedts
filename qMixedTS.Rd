\name{qMixedTS}
\alias{qMixedTS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{MixedTS Quantile function.
%%  ~~function to do ... ~~
}
\description{ Quantile function of mixed Tempered Stable distribution 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
qMixedTS(p, mu0, mu, sig, a, alpha, lambda_p, lambda_m,
         setSup = NULL, setInf = NULL, MixingDens = "Gamma", interval = c(-1000, 1000))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{ vector of probabilities.
%%     ~~Describe \code{p} here~~
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
  \item{setSup}{ Internal parameter.
%%     ~~Describe \code{setSup} here~~
}
  \item{setInf}{ Internal parameter.
%%     ~~Describe \code{setInf} here~~
}
  \item{MixingDens}{ Internal parameter.
%%     ~~Describe \code{MixingDens} here~~
}
  \item{interval}{ Internal parameter.
%%     ~~Describe \code{interval} here~~
}
}
\details{ This function returns the quantiles of mixed tempered stable random variable.
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

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{

x<-seq(-3,3,length=400) 
# Plot cumulative distribution function
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
