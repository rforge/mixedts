\name{setMixedTS.param}
\alias{setMixedTS.param}
\alias{MixedTS-parameters}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Mixed Tempered Stable distribution
%%  ~~function to do ... ~~
}
\description{ \code{setMixedTS} describes the Mixed Tempered Stable distribution introduced in Rroji and Mercuri (2015):

\bold{Definition}

We say that a continuous random variable Y follows a Mixed Tempered Stable distribution if:

    \code{Y= mu0+ mu*V + sigma*sqrt{V}*Z}

The conditional distribution of random variable given V=v is a standardized Tempered Stable with parameters \code{(alpha, lambda_p*sqrt{v}, lambda_m*sqrt{v})} (see Kuchler, U. and Tappe, S. 2014). The distribution of V is infinitely divisible defined on the positive axis.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
setMixedTS.param(mu0 = numeric(), mu = numeric(),
  sigma = numeric(), a, alpha = numeric(),
  lambda_p = numeric(), lambda_m = numeric(),
  param = numeric(), Mixing = "Gamma", paramMixing = list(), Parametrization = "A")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mu0}{ a numeric object. \code{mu0} parameter belongs to the real axis.
}
  \item{mu}{ a numeric object. \code{mu} parameter belongs to the real axis
%%     ~~Describe \code{mu} here~~
}
  \item{sigma}{ a numeric object. \code{sigma} parameter assumes value from zero to infinity.
%%     ~~Describe \code{sigma} here~~
}
  \item{a}{ a vector object. If numeric, the mixing density \code{V} is a Gamma and \code{a} is the value of the shape parameter. If string, \code{a} is the log of the moment generating function of the mixing density \code{V}.
}
  \item{alpha}{ a numeric object that takes value from 0 to 2. If alpha is fixed to 2, the Mixed Tempered Stable becomes the Normal Variance Mean mixture.
%%     ~~Describe \code{alpha} here~~
}
  \item{lambda_p}{ a positive numeric object. It is the right tempering parameter of the random variable \code{X}.
%%     ~~Describe \code{lambda_p} here~~
}
  \item{lambda_m}{ a positive numeric object. It is the left tempering parameter of the random variable \code{X}
%%     ~~Describe \code{lambda_m} here~~
}
  \item{param}{ a numeric object containing the Mixed Tempered Stable parameters. It is not necessary if we use the previous inputs for defining the distribution. See documentation for more details.
%%     ~~Describe \code{param} here~~
}
  \item{Mixing}{ a string object indicating the nature of the mixing density \code{V}. If \code{Mixing="Gamma"} (default value), the \code{V} randm variable is a Gamma. If \code{Mixing="Gamma"}, the user have to specify the log of the moment generating function of the \code{V} random variable.
%%     ~~Describe \code{Mixing} here~~
}
  \item{paramMixing}{ a list object. It is an empty list when \code{Mixing="Gamma"}. If \code{Mixing="User"}, it is used to pass the values of the Mixing density parameters defined by the User through slot \code{a}.
%%     ~~Describe \code{paramMixing} here~~
}
  \item{Parametrization}{ a character string. If \code{Parametrization="A"} the default, we use the following definition for MixedTS with gamma density

      \code{Y= mu0+ mu*V + sqrt{V}*Z}

where \code{V} is distributed as a \code{Gamma(a, sigma^2)}. Otherwise if \code{Parametrization="B"} we have:

      \code{Y= mu0+ mu*V + sigma*sqrt{V}*Z}
  }

where \code{V} is distributed as a \code{Gamma(a, 1)}.

}
\details{
For particular choices of the tempering parameters the tails of the MixedTS distribution can be heavy or semi-heavy. In particular if the Mixing density is a Gamma, we get the Variance Gamma (Madan and Seneta 1990) and the symmetric Geo-Stable distribution as special cases.
}
\value{ This function returns an object of class \code{"\link{param.MixedTS}"}.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
 Barndorff-Nielsen,O.E., Kent,J. and Sorensen, M. (1982): Normal variance-mean mixtures and z-distributions, \emph{International Statistical Review}, 50, 145-159.

Kuchler, U. and Tappe, S. (2014): Exponential stockmodels driven by tempered stable processes. \emph{Journal of Econometrics},181 (1), 53-63.

Madan, D.B. and Seneta E. (1990): The variance gamma (V.G.) model for share market returns, \emph{Journal of Business}, 63, 511-524

Rroji, E and Mercuri, L.(2015): Mixed Tempered Stable distribution \emph{Quantitative Finance}, 15,(9) pages 1559-1569.
}
\note{
This class of distributions has the Normal Variance Mean Mixture (Barndorff-Nielsen et al. 1982) as special case.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
# Mixed Tempered Stable with Gamma Mixing density.

ParamEx1<-setMixedTS.param(mu0=0, mu=0, sigma=0.4, a=1.5,
                           alpha=0.8, lambda_p=4, lambda_m=1)


# Mixed Tempered Stable with Inverse Gaussian Mixing density.
## As first step we set the "a" parameter
## equal to the log mgf of the inverse gaussian random variable
# The log mgf of an Ig with parameter (lamb, mu1) is defined as:

logmgf<-("lamb/mu1*(1-sqrt(1-2*mu1^2/lamb*u))")
Mix<-"User"

# The parameters of the mixing density are set by the following command
# line:

parMix<-list(lamb=1,mu1=1)

ParamEx2<-setMixedTS.param(mu0=0, mu=0, sigma=0.4, a=logmgf,
                           alpha=0.8, lambda_p=4, lambda_m=1,
                           Mixing=Mix,paramMixing=parMix)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
