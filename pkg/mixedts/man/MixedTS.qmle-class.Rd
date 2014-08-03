\name{MixedTS.qmle-class}
\Rdversion{1.1}
\docType{class}
\alias{MixedTS.qmle-class}
\alias{initialize,MixedTS.qmle-method}

\title{\code{MixedTS.qmle}: a class for Maximum Likelihood  of Mixed Tempered Stable}
\description{
This class is constructed by function \code{MixedTS.qmle}. It is a subclass for the \code{MixedTS-class}
}
\section{Objects from the Class}{
Objects can be created by function \code{MixedTS.qmle}.
}

\section{Slots}{
  \describe{
    \item{\code{time}:}{Object of class \code{"numeric"}. Computational Time.}
    \item{\code{coef}:}{Object of class \code{"numeric"}. Estimated parameters.}
    \item{\code{vcov}:}{Object of class \code{"matrix"}. Approximate variance-covariance matrix.}
    \item{\code{min}:}{Object of class \code{"numeric"}.  Minimum value of objective function. }
    \item{\code{details}:}{Object of class \code{"list"}. A list as returned from \code{constrOptim} }
    \item{\code{nobs}:}{Object of class \code{"integer"}. Number of observation. }
    \item{\code{method}:}{Object of class \code{"character"}. The optimization method used.}
    \item{\code{Data}:}{Object of class \code{"numeric"}. See \code{MixedTS-class}.
    }
    \item{\code{dens}:}{Object of class \code{"numeric"}. See \code{MixedTS-class}. }
    \item{\code{prob}:}{Object of class \code{"numeric"}. See \code{MixedTS-class}. }
    \item{\code{xMixedTS}:}{Object of class \code{"numeric"} . See \code{MixedTS-class}. }
    \item{\code{quantile}:}{Object of class \code{"logical"}. See \code{MixedTS-class}. }
    \item{\code{mu0}:}{Object of class \code{"numeric"}. See \code{MixedTS-class}. }
    \item{\code{mu}:}{Object of class \code{"numeric"}. See \code{MixedTS-class}. }
    \item{\code{sigma}:}{Object of class \code{"numeric"}. See \code{MixedTS-class}. }
    \item{\code{a}:}{Object of class \code{"vector"}. See \code{MixedTS-class}. }
    \item{\code{alpha}:}{Object of class \code{"numeric"}. See \code{MixedTS-class}. }
    \item{\code{lambda_p}:}{Object of class \code{"numeric"}. See \code{MixedTS-class}. }
    \item{\code{lambda_m}:}{Object of class \code{"numeric"}. See \code{MixedTS-class}. }
    \item{\code{Mixing}:}{Object of class \code{"character"}. See \code{MixedTS-class}. }
    \item{\code{paramMixing}:}{Object of class \code{"list"} . See \code{MixedTS-class}.}
    \item{\code{MixingLogMGF}:}{Object of class \code{"OptionalFunction"}. See \code{MixedTS-class}. }
  }
}
\section{Extends}{
Class \code{"\linkS4class{MixedTS}"}, directly.
Class \code{"\linkS4class{param.MixedTS}"}, by class "MixedTS", distance 2.
}
\section{Methods}{
  \describe{
    \item{summary}{\code{signature(.Object = "MixedTS.qmle")}}
    \item{coef}{\code{signature(.Object = "MixedTS.qmle")}}
    \item{vcov}{\code{signature(.Object = "MixedTS.qmle")}}
    \item{logLik}{\code{signature(.Object = "MixedTS.qmle")}}
    \item{BIC}{\code{signature(.Object = "MixedTS.qmle")}}
    \item{AIC}{\code{signature(.Object = "MixedTS.qmle")}}
   }
}
