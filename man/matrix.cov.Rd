\name{matrix.cov}
\alias{matrix.cov}


\title{Covariance}
\description{For a single gene, computes the transformed or untransformed sample covariance matrix 
if one biological condition, or pooled sample covariance matrix if two or more biological conditions.}
\usage{
matrix.cov(x, k, trans = TRUE, c.grp = NULL, use = "complete.obs")
}
\arguments{
   \item{x}{a numeric vector giving the log-ratios or log-values for a gene, sorted 
           in ascending order by biological condition, replicate, and time groups.}
  \item{k}{a positive integer giving the number of time points.}
  \item{trans}{logical. Should the Helmert transformation be performed?}
  \item{c.grp}{a numeric vector corresponding to the biological condition group for each element of   
\code{x}.}
  \item{use}{character. The same as the \code{use} in stats function \code{cov}. 
             The default uses complete observations.}
}

\details{
This function is for internal use only and is not to be called by the user.
}

\value{
A numeric matrix. 
}

\references{
Becker, R. A., Chambers, J. M. and Wilks, A. R. (1988) \emph{The New S
     Language}. Wadsworth & Brooks/Cole.
}

\seealso{
\code{\link{cov}}, \code{\link{ot.helmert}}.
}

\author{Yu Chuan Tai \email{yuchuan@stat.berkeley.edu}}
\keyword{internal}

