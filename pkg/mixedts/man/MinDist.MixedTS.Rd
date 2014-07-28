\name{MinDist.MixedTS}
\alias{MinDist.MixedTS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Fitting Mixed Tempered Stable Density
}
\description{The function returns the parameters that minimize the mean squared error between empirical and theoretical density.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MinDist.MixedTS(data, param0, method = "L-BFGS-B", fixed.param = NULL, 
        lower.param = NULL, upper.param = NULL, MixingDens = "Gamma", N = 50)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ ...
%%     ~~Describe \code{data} here~~
}
  \item{param0}{ ...
%%     ~~Describe \code{param0} here~~
}
  \item{method}{ ...
%%     ~~Describe \code{method} here~~
}
  \item{fixed.param}{ ...
%%     ~~Describe \code{fixed.param} here~~
}
  \item{lower.param}{ ...
%%     ~~Describe \code{lower.param} here~~
}
  \item{upper.param}{ ...
%%     ~~Describe \code{upper.param} here~~
}
  \item{MixingDens}{ ...
%%     ~~Describe \code{MixingDens} here~~
}
  \item{N}{ Number of Class.
%%     ~~Describe \code{N} here~~
}
}

\value{ Estimated parameters and the objective function
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
