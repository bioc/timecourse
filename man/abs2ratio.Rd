\name{abs2ratio}
\alias{abs2ratio}

\title{Convert log-values to log-ratios}
\description{For a single gene, computes the log ratios between time courses 
from two paired biological conditions.}
\usage{
abs2ratio(x, mn, k, c.grp, reference)
}
\arguments{
  \item{x}{a numeric vector giving the log-values of a gene with two paired biological conditions,
           sorted in ascending order by biological condition, replicate, and time groups.}
  \item{mn}{a numeric matrix giving the sample sizes for the two biological conditions.}
  \item{k}{a positive integer giving the number of time points.}
  \item{c.grp}{an numeric or character vector with length equals to that of \code{x},
                 giving the biological condition group for each element of \code{x}.}
  \item{reference}{a numeric value or character assigning the reference biological condition.}
}
\value{
a numeric vector containing log-ratios between two paired biological conditions.
}

\details{
This function is for internal use only and is not to be called by the user.

}
\seealso{
\code{\link{mb.paired}}.
}

\author{Yu Chuan Tai \email{yuchuan@stat.berkeley.edu}}
\keyword{internal}

