\title{Gene Temporal Profile Plot}
\name{plotProfile}
\alias{plotProfile}
\description{
Plots the longitudinal temporal profile of a gene.
}
\usage{
plotProfile(object, stats=c("HotellingT2", "MB"), ranking=1, gid=NULL, gnames=NULL, desc=NULL,
type=c("p","l","b"), col=2:100, lty=1:100, pch=1:100, lwd=2, xlab="Time",
ylab="Expression", legloc=NULL, xlim=NULL, ylim=NULL, cex.main=1,...)
}
\arguments{
  \item{object}{a \code{MArrayTC} object.}
  \item{stats}{a character indicating which statistic the \code{ranking} is based on.}
  \item{ranking}{a numeric value giving the ranking of the gene to be plotted.}
  \item{gid}{an optional character giving the ID of the gene to be plotted.}
  \item{gnames}{an optional character vector with the \eqn{i_th} element corresponds to the gene ID of the \eqn{i_th} gene in \code{object$M}.}
  \item{desc}{an optional character vector with the \eqn{i_th} element corresponds to the gene description of the \eqn{i_th} gene in \code{object$M}.}
  \item{type}{a character indicating the plot type, \code{"p"} for points, \code{"l"} for lines, and \code{"b"} for both.}
  \item{col}{a character or numeric vector giving the colors for different biological conditions. Default is \code{2:100}.}
  \item{lty}{a character or numeric vector giving the line types for different replicates. Default is \code{1:100}.}
  \item{pch}{a character or numeric vector giving the point types for different replicates. Default is \code{1:100}.}
  \item{lwd}{optional. The default sets to 2.}
  \item{xlab}{character. The label for the x-axis.}
  \item{ylab}{character. The label for the y-axis.} 
  \item{legloc}{an optional vector giving the location of the legend.}
  \item{xlim}{an optional vector giving the upper- and lower- limits of x-axis.}
  \item{ylim}{an optional vector giving the upper- and lower- limits of y-axis.}
  \item{cex.main}{optional. The default sets to 1}
  \item{...}{any other arguments passed onto \code{plot}}
}

\details{
This function takes an object of \code{MArrayTC} as the input and plots
the temporal profile of a single gene. The user can specify either the
ranking based on \code{stats} or the gene ID of the gene to be plotted.   

See \code{\link[graphics]{points}} for possible values for \code{pch},
\code{col} and \code{cex}. 

See \code{mb.long} for examples.
}


\author{Yu Chuan Tai \email{yuchuan@stat.berkeley.edu}}



\keyword{hplot}

