\name{mb.MANOVA}
\alias{mb.MANOVA}

\title{Multivariate Empirical Bayes Analysis of Variance for Longitudinal Replicated 
Developmental Microarray Time Course Data}
\description{Computes the MB-statistics for longitudinal replicated developmental microarray time course 
data with multiple biological conditions.}
\usage{
mb.MANOVA(object, times, D, size, nu = NULL, Lambda = NULL, 
beta.d = NULL, beta = NULL, alpha.d = NULL, alpha = NULL, 
condition.grp, time.grp = NULL, rep.grp = NULL, p = 0.02)
}
\arguments{
  \item{object}{Required. An object of class \code{matrix}, \code{MAList},
    \code{marrayNorm}, or \code{exprSet} containing log-ratios or
    log-values of expression for a series of microarrays.} 
  \item{times}{Required. A positive integer giving the number of time points.}
  \item{D}{Required. A positive integer giving the number of biological conditions. \code{D>1}}
  \item{size}{Required. A numeric matrix corresponding to the sample sizes for all genes across
          different biological conditions, when biological conditions are sorted in ascending order.
          Rows represent genes while columns represent biological conditions.}
  \item{nu}{an optional positive value giving the degrees of moderation for the fully moderated Wilks' 
Lambda.}  
  \item{Lambda}{an optional numeric matrix giving the common covariance matrix for the fully
   moderated Wilks' Lambda.}
  \item{beta.d}{an optional numeric vector of length \code{D} giving the condition-specific scale 
parameters for the common covariance matrix of the expected time course vectors under the alternative.}
  \item{beta}{an optional numeric value giving the scale parameter for the common covariance matrix of the common expected time course vector under the null.}
  \item{alpha.d}{an optional numeric matrix giving the condition-specific means of the expected
  time course vectors under the alternative.}
  \item{alpha}{an optional numeric vector of length \code{times} giving the common mean of the 
 expected time course vector under the null.}
\item{condition.grp}{Required. A numeric or character vector with length equals to the number of arrays,
assigning the biological condition group to each array.}
  \item{rep.grp}{an optional numeric or character vector with length equals to the number of arrays,
assigning the replicate group to each array.}
  \item{time.grp}{an optional numeric vector with length equals to the number of arrays,
assigning the time point group to each array.}
  \item{p}{a numeric value between 0 and 1, assumed proportion of genes
which are differentially expressed.}
}
\value{

Object of \code{MArrayTC}.
  
}
\details{
  This function implements the multivariate empirical Bayes Analysis
  of Variance model for identifying genes with different temporal
  profiles across multiple biological conditions,
  as described in Tai (2005).  
  
}

\seealso{
timecourse Vignette.
}
\author{Yu Chuan Tai  \email{yuchuan@stat.berkeley.edu}}
\references{
Yu Chuan Tai (2005). Multivariate empirical Bayes models for replicated microarray time course data.
Ph.D. dissertation. Division of Biostatistics, University of California, Berkeley.

Yu Chuan Tai and Terence P. Speed (2005). Statistical analysis of microarray time course data. In: DNA 
Microarrays, U. Nuber (ed.), BIOS Scientific Publishers Limited, Taylor & Francis, 4 Park Square, Milton 
Park, Abingdon OX14 4RN, Chapter 20. 
}
\examples{


SS <- matrix(c(    0.01, -0.0008,   -0.003,     0.007,  0.002,
                -0.0008,    0.02,    0.002,   -0.0004, -0.001,
                 -0.003,   0.002,     0.03,   -0.0054, -0.009,
                  0.007, -0.0004, -0.00538,      0.02, 0.0008,
                  0.002,  -0.001,   -0.009,    0.0008,  0.07), ncol=5)


sim.Sigma <- function()
{
   S <- matrix(rep(0,25),ncol=5)
   x <- mvrnorm(n=10, mu=rep(0,5), Sigma=10*SS)
   for(i in 1:10)
       S <- S+crossprod(t(x[i,]))

   solve(S)

}

## Now let's simulate a dataset with three biological conditions
## 500 genes in total, 10 of them have different expected time course profiles
## across biological conditions
## the first condition has 3 replicates, while the second condition has 4 replicates, 
## and the third condition has 2 replicates. 5 time points for each condition.

sim.data <- function(x, indx=1)
{
   mu <- rep(runif(1,8,x[1]),5)
   if(indx==1)
     res <- c(as.numeric(t(mvrnorm(n=3, mu=mu+rnorm(5,sd=5), Sigma=sim.Sigma()))),
             as.numeric(t(mvrnorm(n=4, mu=mu+rnorm(5,sd=3.2), Sigma=sim.Sigma()))),
             as.numeric(t(mvrnorm(n=2, mu=mu+rnorm(5,sd=2), Sigma=sim.Sigma()))))

   if(indx==0) res <- as.numeric(t(mvrnorm(n=9, mu=mu+rnorm(5,sd=3), Sigma=sim.Sigma())))
   res
}

M <- matrix(rep(14,500*45), ncol=45)
M[1:10,] <- t(apply(M[1:10,],1,sim.data))
M[11:500,] <- t(apply(M[11:500,],1,sim.data, 0))


assay <- rep(c("1.2.04","2.4.04","3.5.04","5.21.04","7.17.04","9.10.04","12.1.04","1.2.05","4.1.05"),each=5)
trt <- c(rep(c("wildtype","mutant1"),each=15),rep("mutant1",5), rep("mutant2", 10))

# Caution: since "mutant1" < "mutant2" < "wildtype", the sample sizes should be in the order of 4,2,3, 
# but NOT 3,4,2. 
size <- matrix(c(4,2,3), byrow=TRUE, nrow=500, ncol=3)
MB.multi <- mb.MANOVA(M, times=5, D=3, size=size, rep.grp=assay, condition.grp=trt)

plotProfile(MB.multi, stats="MB", type="b") # plots the no. 1 gene
}
\keyword{multivariate}
