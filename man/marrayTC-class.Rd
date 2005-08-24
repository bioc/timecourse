\name{MArrayTC-class}
\docType{class}
\alias{MArrayTC-class}
\alias{LargeDataObject-class}
\title{Microarray Time Course Object- class}

\description{
A list-based class for storing the analysis results from the multivariate
empirical Bayes models of differential expression for longitudinal replicated
developmental microarray time course data.
Objects are normally created by \code{\link{mb.long}} and \code{\link{mb.MANOVA}}.
}

\section{Slots/Components}{
  \code{MArrayTC} objects do not contain any slots (apart from
  \code{.Data}) but they should contain the following list components:
  \describe{
    \item{\code{M}:}{input \code{matrix} of log-ratios or log-values of
      expression for a series of microarrays.}
  }
  Objects may also contain the following optional components:
  \describe{
     \item{\code{prop}:}{\code{numeric} value giving the proportion of
      differentially expressed genes.}
    \item{\code{nu}:}{\code{numeric} value containing the estimated
      amount of moderation.}
    \item{\code{Lambda}:}{the estimated Lambda.}
    \item{\code{Lambda1}:}{the estimated Lambda1.}
    \item{\code{eta}:}{the estimated prior scale parameter.}
    \item{\code{alpha}:}{the estimated common
          mean of the expected time course vector under the null.}
     \item{\code{alpha.d}:}{the estimated condition-specific     
          means of the expected time course vectors under the alternative.}
     \item{\code{beta}:}{the estimated scale parameter for the
          common covariance matrix of the common expected time course
          vector under the null.}
     \item{\code{beta.d}:}{the estimated condition-specific scale parameters for the common
          covariance matrix of the expected time course vectors under
          the alternative.} 
    \item{\code{percent}:}{\code{numeric} matrix containing the percent of moderation 
     corresponding to each sample size for the longitudinal one- and two- sample problems.}
    \item{\code{size}:}{\code{numeric} vector or matrix containing the sample sizes for 
    all genes corresponding to different biological conditions, when the latter are sorted in ascending 
    order.}
    \item{\code{con.group}:}{\code{numeric} or \code{character} vector giving the biological
    condition group of each array. The \eqn{i_th} element of \code{con.group} corresponds to the
    biological condition of the \eqn{i_th} column of \code{M}.}
    \item{\code{rep.group}:}{\code{numeric} or \code{character} vector giving the replicate
    group of each array. The \eqn{i_th} element of \code{rep.group} corresponds to the 
    replicate of the \eqn{i_th} column of \code{M}.}
    \item{\code{time.group}:}{\code{numeric} vector
      giving the time group of each array. The \eqn{i_th} element of \code{time.group} corresponds to the 
    time of the \eqn{i_th} column of \code{M}.}
    \item{\code{HotellingT2}:}{\code{numeric} vector giving the 
    \eqn{\tilde{T}^2} statistics of differential expression.}
    \item{\code{MB}:}{\code{numeric} vector giving the MB-statistics of differential expression.}
    \item{\code{pos.HotellingT2}:}{\code{numeric} vector whose \eqn{i_th} element corresponds to
    the index of the gene with ranking \eqn{i} in \code{HotellingT2}.}
    \item{\code{pos.MB}:}{\code{numeric} vector whose \eqn{i_th} element corresponds to
    the index of the gene with ranking \eqn{i} in \code{MB}.}
    \item{\code{geneNames}:}{\code{character} vector giving gene names.}
    \item{\code{descriptions}:}{\code{character} vector giving gene descriptions.}
  }
}

\section{Methods}{ 
MArrayTC extends the 
\code{\link[limma:LargeDataObject]{LargeDataObject}} class in package limma, and
inherits a \code{show} method from there.

The function \code{plotProfile} takes \code{MArrayTC} as the input argument.
}

\author{Yu Chuan Tai \email{yuchuan@stat.berkeley.edu}}


\keyword{classes}

