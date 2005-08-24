\name{univ.func}
\alias{univ.func}

\title{Univariate Data}
\description{Transforms multivariate vectors into univariate values using the Helmert matrix.}
\usage{
univ.func(dummy, M, k, n, indx = 1)
}
\arguments{
  \item{dummy}{a numeric gene index.}
  \item{M}{a numeric matrix containing the log-values or log-ratios of a gene.}
  \item{k}{a positive integer giving the number of time points.}
  \item{n}{a positive integer giving the number of replicates.}
  \item{indx}{a positive integer between 1 and k, indicating which row of the Helmert
  matrix to transform the vectors.}
}
\value{
A numeric vector with length equals to \code{n}.
}

\details{
This function is for internal use only and is not to be called by the user.
}

\seealso{
\code{\link{ot.helmert}}.
}


\author{Yu Chuan Tai \email{yuchuan@stat.berkeley.edu}}
\keyword{internal}

