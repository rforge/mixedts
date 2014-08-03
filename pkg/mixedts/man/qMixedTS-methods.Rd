\name{qMixedTS-methods}
\docType{methods}
\alias{qMixedTS}
\alias{qMixedTS-methods}
\alias{qMixedTS,param.MixedTS-method}
\title{Quantile of Mixed Tempered Stable distribution}
\description{
 This Method returns the quantile of a Mixed Tempered Stable.
}
\section{Methods}{
\describe{

\item{\code{signature(object = "param.MixedTS",x = numeric(), setSup=NULL,setInf=NULL,N=2^10)}}{
This method returns an object of class \code{MixedTS} where the slot \code{prob} contains the value of the quantile evaluated on the \code{x} (x is the probability). \code{setSup} and \code{setInf} are used to choose \code{+ infinity} and \code{- infinty}. \code{N} is the number of point used for discretization in \code{fft} algorithm.
}
}}
