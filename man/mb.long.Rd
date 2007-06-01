\name{mb.long}
\alias{mb.long}
\alias{mb.1D}
\alias{mb.paired}
\alias{mb.2D}

\title{Multivariate Empirical Bayes Statistics for Longitudinal Replicated Developmental Microarray Time 
Course Data}
\description{Computes the \eqn{\tilde{T}^2} statistics and/or the 
MB-statistics of differential expression for longitudinal replicated developmental microarray time course 
data by multivariate empirical Bayes shrinkage of gene-specific sample variance-covariance matrices
towards a common matrix.}
\usage{
mb.long(object, method = c("1D", "paired", "2D"), type = c("none", "robust"), 
times, reps, prior.df = NULL, prior.COV = NULL, 
prior.eta = NULL, condition.grp = NULL, rep.grp = NULL, time.grp = NULL, 
one.sample = FALSE, ref = NULL, p = 0.02, out.t = FALSE, 
tuning = 1.345, HotellingT2.only=TRUE)
}
\arguments{
  \item{object}{Required. An object of class \code{matrix}, \code{MAList},
    \code{marrayNorm}, or \code{exprSet} containing log-ratios or
    log-values of expression for a series of microarrays.}
  \item{method}{a character string, \code{"1D"} for the one-sample case where genes
                of interest are those which change over time,  
                \code{"paired"} for the one-sample case where genes of interest are
		those whose expected temporal profiles do not stay 0, for example, cDNA microarrays, 
                or the paired two-sample case where genes of interest are those with different
                expected temporal profiles across 2 biological conditions,
                \code{"2D"} for the independent two-sample case where genes of interest are those
                with different expected temporal profiles across 2 biological conditions.
                The default is \code{"1D"}.}
  \item{type}{a character string, indicating whether possible outliers should be down-weighted.} 
  \item{times}{Required. A positive integer giving the number of time points.}
  \item{reps}{Required. A numeric vector or matrix corresponding to the sample sizes for all 
          genes across different biological conditions, when biological conditions are sorted in 
          ascending order. If a matrix, rows represent genes while columns represent biological conditions. }
  \item{prior.df}{an optional positive value giving the degrees of moderation.}  
  \item{prior.COV}{an optional numeric matrix giving the common covariance matrix
    to which the gene-specific sample covariances are smoothed toward.}  
  \item{prior.eta}{an optional numeric value giving the scale parameter for the covariance matrix for
   the expected time course profile.}
  \item{condition.grp}{a numeric or character vector with length equals to the number of arrays,
                       assigning the biological condition group of each array. Required if 
                       \code{method=2D}.}
  \item{rep.grp}{an optional numeric or character vector with length equals to the number of arrays,
                 assigning the replicate group of each array.}
  \item{time.grp}{an optional numeric vector with length equals to the number of arrays,
                 assigning the time point group of each array.}
  \item{one.sample}{Is it a one-sample problem? Only specify this argument when \code{method=paired}. 
                    The default is \code{FALSE} which means it is a paired two-sample problem.}
  \item{ref}{an optional numeric value or character specifying the name
          of reference biological condition. The default uses the
          first element of \code{condition.grp}. Only specify this argument when
          \code{method=paired} and \code{one.sample} is \code{FALSE}.}
  \item{p}{a numeric value between 0 and 1, assumed proportion of genes
          which are differentially expressed.}
  \item{out.t}{logical. Should the moderated multivariate t-statistics be outputed? The default is 
  \code{FALSE}.}
  \item{tuning}{the tuning constant for the Huber weight function with a default 1.345.}
  \item{HotellingT2.only}{logical. Should only the HotellingT2 statistics be outputed? This should be 
set as \code{TRUE} (default) when the sample size(s) are the same across genes, in order to reduce
computational time. }
}
\value{

Object of \code{MArrayTC}.
  
}
\details{
  This function implements the multivariate empirical Bayes statistics
  described in Tai and Speed (2004), to rank genes in the order of
  interest from longitudinal replicated developmental microarray time course
  experiments. It calls one of the following functions, 
  depending
  on which \code{method} is used: \code{mb.1D}, 
  \code{mb.paired}, and \code{mb.2D}.

  
  The arguments \code{condition.grp}, \code{rep.grp}, and
  \code{time.grp}, if specified, should have lengths equal to the number
  of arrays. The \eqn{i_th} elements of these three arguments should
  correspond to the biological condition, replicate, and time for the \eqn{i_th} column (array) in the
  expression value matrix of the input object, respectively.
  The default assumes the columns of \code{M} are in the ascending order of \code{condition.grp} first,
  and then \code{rep.grp}, and finally \code{time.grp}. 
 
  Arguments \code{one.sample} and \code{ref} are for \code{method=paired} only.

  When \code{type=robust}, the numerator of the \eqn{\tilde{T}^2} statistic is calculated using
  the weighted average time course vector(s), where the weight at each data point 
  is determined using Huber's weight function with the default tuning constant 1.345.

  Warning: When there are only 2 replicates within conditions,
           \code{type="robust"} produces the same rankings as \code{type="none"} 
           since there is no concensus on gene expression values. 
           Check the output weights for these outliers.   
}

\seealso{
timecourse Vignette.
}
\author{Yu Chuan Tai  \email{yuchuan@stat.berkeley.edu}}
\references{
Yu Chuan Tai and Terence P. Speed (2006). A multivariate empirical Bayes statistic
for replicated microarray time course data. Annals of Statistics 34(5):2387-2412. 

Yu Chuan Tai and Terence P. Speed (2005). Statistical analysis of microarray time course data. In: DNA 
Microarrays, U. Nuber (ed.), BIOS Scientific Publishers Limited, Taylor & Francis, 4 Park Square, Milton 
Park, Abingdon OX14 4RN, Chapter 20. 


P. J. Huber (2004). Robust Statistics. 
\emph{Wiley series in probability and mathematical statistics}.
}
\examples{
data(fruitfly)
colnames(fruitfly) ## check if arrays are arranged in the default order
gnames <- rownames(fruitfly)
assay <- rep(c("A", "B", "C"), each = 12)
time.grp <- rep(c(1:12), 3)
size <- rep(3, nrow(fruitfly))

out1 <- mb.long(fruitfly, times=12, reps=size, rep.grp = assay, time.grp = time.grp)
summary(out1)
plotProfile(out1, type="b", gnames=gnames, legloc=c(2,15), pch=c("A","B","C"), xlab="Hour")

## Simulate gene expression data
## Note: this simulation is for demonstration purpose only,
## and does not necessarily reflect the real 
## features of longitudinal time course data

## one biological condition, 5 time points, 3 replicates
## 500 genes, 10 genes change over time


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

sim.data1 <- function(x, indx=1)
{
   mu <- rep(runif(1,8,x[1]),5)
   if(indx==1) res <- as.numeric(t(mvrnorm(n=3, mu=mu+rnorm(5,sd=4), Sigma=sim.Sigma())))
   if(indx==0) res <- as.numeric(t(mvrnorm(n=3, mu=mu, Sigma=sim.Sigma())))
   res
}

M1 <- matrix(rep(14,500*15), ncol=15)
M1[1:10,] <- t(apply(M1[1:10,],1,sim.data1))
M1[11:500,] <- t(apply(M1[11:500,],1,sim.data1, 0))

## Which genes are nonconstant?
MB.1D1 <- mb.long(M1, times=5, reps=rep(3, 500))
MB.1D1$percent  # check the percent of moderation

plotProfile(MB.1D1,type="b") # plots the no. 1 gene
plotProfile(MB.1D1,type="b",ranking=10) # plots the no. 10 gene
genenames <- as.character(1:500)
plotProfile(MB.1D1, type="b", gid="8", gnames=genenames) #plots the gene with ID "8"

## 
MB.1D1.r <- mb.long(M1, type="r", times=5, reps=rep(3, 500))
plotProfile(MB.1D1.r,type="b",gnames=genenames)
plotProfile(MB.1D1.r,type="b", gid="1", gnames=genenames) #plots the gene with ID "1" 

## assign the following labellings to columns of M1
## which is actually the same as the default
## Not Run
trt <- rep("wildtype", 15)
assay <- rep(c("A","B","C"), rep(5,3))
time.grp <- rep(c(0, 1, 3, 4, 6), 3)

## MB.1D2 should give the same results as MB.1D1
#MB.1D2 <- mb.long(M1, times=5, reps=rep(3, 500), condition.grp = trt, rep.grp = assay, 
#time.grp=time.grp)

## suppose now the replicates are in this order instead
assay <- rep(c("A","C","B"), rep(5,3))

## then
MB.1D3 <- mb.long(M1, times=5, reps=rep(3, 500), condition.grp = trt, rep.grp = assay, time.grp=time.grp)
MB.1D3$rep.group  #check the replicate and time group
MB.1D3$time.group


## Now let's simulate another dataset with two biological conditions
## 500 genes also, 10 of them have different expected time course profiles
## between these two biological conditions  
## 3 replicates, 5 time points for each condition

sim.data2 <- function(x, indx=1)
{
   mu <- rep(runif(1,8,x[1]),5)
   if(indx==1)
     res <- c(as.numeric(t(mvrnorm(n=3, mu=mu+rnorm(5,sd=5), Sigma=sim.Sigma()))),
             as.numeric(t(mvrnorm(n=3, mu=mu+rnorm(5,sd=3.2), Sigma=sim.Sigma()))))

   if(indx==0) res <- as.numeric(t(mvrnorm(n=6, mu=mu+rnorm(5,sd=3), Sigma=sim.Sigma())))
   res 
}

M2 <- matrix(rep(14,500*30), ncol=30)
M2[1:10,] <- t(apply(M2[1:10,],1,sim.data2))
M2[11:500,] <- t(apply(M2[11:500,],1,sim.data2, 0))

## assume it is a paired two-sample problem
trt <- rep(c("wt","mt"),each=15)
assay <- rep(rep(c("1.2.04","2.4.04","3.5.04"),each=5),2)
size <- matrix(3, nrow=500, ncol=2)
MB.paired <- mb.long(M2, method="paired", times=5, reps=size, condition.grp=trt, rep.grp=assay)
MB.paired$con.group # check the condition, replicate and time groups
MB.paired$rep.group
MB.paired$time.group

plotProfile(MB.paired, type="b")
genenames <- as.character(1:500)
plotProfile(MB.paired, gid="12", type="b", gnames=genenames) #plots the gene with ID "12"

### assume it is a unpaired two-sample problem
assay <- rep(c("1.2.04","2.4.04","3.5.04","5.21.04","7.17.04","8.4.04"),each=5)
MB.2D <- mb.long(M2, method="2", times=5, reps=size, condition.grp=trt, rep.grp=assay)
MB.2D$con.group # check the condition, replicate and time groups
MB.2D$rep.group
MB.2D$time.group 

plotProfile(MB.2D,type="b", gnames=genenames) # plot the no. 1 gene


## Now let's simulate another dataset with two biological conditions
## 500 genes also, 10 of them have different expected time course profiles
## between these two biological conditions
## the first condition has 3 replicates, while the second condition has 4 replicates, 
## 5 time points for each condition

sim.data3 <- function(x, indx=1)
{
   mu <- rep(runif(1,8,x[1]),5)
   if(indx==1)
     res <- c(as.numeric(t(mvrnorm(n=3, mu=mu+rnorm(5,sd=5), Sigma=sim.Sigma()))),
             as.numeric(t(mvrnorm(n=4, mu=mu+rnorm(5,sd=3.2), Sigma=sim.Sigma()))))

   if(indx==0) res <- as.numeric(t(mvrnorm(n=7, mu=mu+rnorm(5,sd=3), Sigma=sim.Sigma())))
   res
}

M3 <- matrix(rep(14,500*35), ncol=35)
M3[1:10,] <- t(apply(M3[1:10,],1,sim.data3))
M3[11:500,] <- t(apply(M3[11:500,],1,sim.data3, 0))

assay <- rep(c("1.2.04","2.4.04","3.5.04","5.21.04","7.17.04","9.10.04","12.1.04"),each=5)
trt <- c(rep(c("wildtype","mutant"),each=15),rep("mutant",5))
## Note that "mutant" < "wildtype", the sample sizes are (4, 3)
size <- matrix(c(4,3), nrow=500, ncol=2, byrow=TRUE)
MB.2D.2 <- mb.long(M3, method="2", times=5, reps=size, rep.grp=assay, condition.grp=trt)
MB.2D.2$con.group # check the condition, replicate and time groups
MB.2D.2$rep.group
MB.2D.2$time.group 

plotProfile(MB.2D.2, type="b") # plot the no. 1 gene
}
\keyword{multivariate}
