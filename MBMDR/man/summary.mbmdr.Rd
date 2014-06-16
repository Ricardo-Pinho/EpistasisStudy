\name{summary.mbmdr}
\Rdversion{2.1}
\alias{summary.mbmdr}
\alias{print.summary.mbmdr}

\title{
  Summary Method for Model Based Multifactor Dimensionality Reduction
}
\description{
  A summary method for mbmdr.
}
\usage{
 \method{summary}{mbmdr}(object,sig.level=0.05,...)
}

\arguments{
  \item{object}{An mbmdr object returned from \code{mbmdr} function.}
  \item{sig.level}{Threshold value for significant p-values.}
  \item{...}{Further arguments passed to or from other methods.}
}

\value{
	It returns a data.frame with those models with a MIN.P value smaller than \code{sig.level}.
}

\seealso{
  \code{\link{mbmdr}}
}
