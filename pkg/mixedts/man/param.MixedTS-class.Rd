\name{param.MixedTS-class}
\Rdversion{1.1}
\docType{class}
\alias{param.MixedTS-class}
\alias{param.MixedTS}
\title{\code{"param.MixedTS"}: A mathematical Description of the Mixed Tempered Stable}
\description{
Main class of the package \code{MixedTS}.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{setMixedTS}.
}
\section{Slots}{
  \describe{
    \item{\code{mu0}:}{a numeric object. mu0 parameter belongs to the real axis. }
    \item{\code{mu}:}{a numeric object. mu parameter belongs to the real axis }
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
  \item{Mixing}{ a string object indicating the nature of the mixing density \code{V}. If \code{Mixing="Gamma"} (default value), the \code{V} randm variable is a Gamma. If \code{Mixing="Gamma"}, the user have to specify the log of the moment generating function of the \code{V} random variable.
%%     ~~Describe \code{Mixing} here~~
}
  \item{paramMixing}{ a list object. It is an empty list when \code{Mixing="Gamma"}. If \code{Mixing="User"}, it is used to pass the values of the Mixing density parameters defined by the User through slot \code{a}.
%%     ~~Describe \code{paramMixing} here~~
}
    \item{\code{MixingLogMGF}:}{This slot contains a function that returns the logarithm of mgf for the Mixing density. The function is  built internally using the information contains into the slots \code{a}, \code{paramMixing}. }
  }
}
\section{Methods}{
  \describe{
    \item{dMixedTS}{\code{signature(object = "param.MixedTS")}: Method for computing density of MixedTS. See \code{"\link{dMixedTS-methods}"} for more details. }
    \item{pMixedTS}{\code{signature(object = "param.MixedTS")}: Method for computing probability of MixedTS. See \code{"\link{pMixedTS-methods}"} for more details. }
    \item{qMixedTS}{\code{signature(object = "param.MixedTS")}: Method for computing quantile of MixedTS. See \code{"\link{qMixedTS-methods}"} for more details. }
    \item{rMixedTS}{\code{signature(object = "param.MixedTS")}: Method for computing random numbers of MixedTS. See \code{"\link{rMixedTS-methods}"} for more details. }
    \item{initialize}{\code{signature(object = "param.MixedTS").}}
	 }
}
