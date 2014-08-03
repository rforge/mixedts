\name{MixedTS-class}
\Rdversion{1.1}
\docType{class}
\alias{MixedTS-class}
\alias{initialize,MixedTS-method}
\alias{plot,MixedTS,ANY-method}

\title{\code{"MixedTS"}: A class for informations about Mixed Tempered Stable}
\description{
Mathematical description of the Mixed Tempered Stable distribution. 

This class inherits from the class \code{param.MixedTS} and is a superclass for \code{MixedTS.qmle-class}.

}
\section{Objects from the Class}{
This object is built by the following methods:

\code{dMixedTS}, \code{pMixedTS}, \code{qMixedTS}, \code{rMixedTS}.
}
\section{Slots}{
  \describe{
    \item{\code{Data}:}{Object of class \code{"numeric"} containing a random number. This slot is filled when the method \code{rMixedTS} is used. }
    \item{\code{dens}:}{Object of class \code{"numeric"} that contains the density of the MixedTS. This slot is filled by \code{dMixedTS}.}
    \item{\code{prob}:}{Object of class \code{"numeric"} that contains the probability of the MixedTS. This slot is filled by \code{pMixedTS} and \code{pMixedTS}. }
    \item{\code{xMixedTS}:}{Object of class \code{"numeric"} that contains the support for the density and probability.}
    \item{\code{quantile}:}{Object of class \code{"logical"}. If \code{TRUE} the object is built by the method \code{qMixedTS}. If \code{FALSE} the object is built by the method \code{qMixedTS}.}
    \item{\code{mu0}:}{Object of class \code{"numeric"}. See \code{param.MixedTS}. }
    \item{\code{mu}:}{Object of class \code{"numeric"}. See \code{param.MixedTS}. }
    \item{\code{sigma}:}{Object of class \code{"numeric"}. See \code{param.MixedTS}.}
    \item{\code{a}:}{Object of class \code{"vector"}. See \code{param.MixedTS}.}
    \item{\code{alpha}:}{Object of class \code{"numeric"}. See \code{param.MixedTS}.}
    \item{\code{lambda_p}:}{Object of class \code{"numeric"}. See \code{param.MixedTS}.}
    \item{\code{lambda_m}:}{Object of class \code{"numeric"}. See \code{param.MixedTS}.}
    \item{\code{Mixing}:}{Object of class \code{"character"}. See \code{param.MixedTS}.}
    \item{\code{paramMixing}:}{Object of class \code{"list"}. See \code{param.MixedTS}.  }
    \item{\code{MixingLogMGF}:}{Object of class \code{"OptionalFunction"}. See \code{param.MixedTS}.}
  }
}
\section{Extends}{
Class \code{"\linkS4class{param.MixedTS}"}, directly.
}
\section{Methods}{
  \describe{
    \item{plot}{\code{signature(x = "MixedTS", ...)}}
	 }
}
