\name{rMixedTS-methods}
\docType{methods}
\alias{rMixedTS}
\alias{rMixedTS-methods}
\alias{rMixedTS,param.MixedTS-method}
\title{ Random number of Mixed Tempered Stable distribution}
\description{
 This Method returns the quantile of a Mixed Tempered Stable.
}
\section{Methods}{
\describe{
\item{\code{signature(object = "param.MixedTS",x = numeric(), setSup=NULL,setInf=NULL,N=2^10)}}{
This method returns an object of class \code{MixedTS} where the slot \code{Data} contains a set of size \code{x} of random numbers. \code{setSup} and \code{setInf} are used to choose \code{+ infinity} and \code{- infinty}. \code{N} is the number of point used for discretization in \code{fft} algorithm.
}
}}
