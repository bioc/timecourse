## (c) 2004-2005 Yu Chuan Tai
## University of California, Berkeley
## Created: August 2004
## Last Modified: Feb. 15, 2005

## Helmert Matrix

ot.helmert <- function(k)
{

    if(missing(k)) stop("The number of time points is missing.")

    if (is.numeric(k) && length(k) == 1)
          if(k > trunc(k)) stop("The number of time points is not an integer.")

    
    levels <- 1:k

    T0 <- matrix(rep(1/sqrt(k), k), byrow=TRUE, ncol=k)

    T1 <- matrix(rep(0,(k-1)*k), ncol=k, byrow=TRUE)
    T1 <- array(1/sqrt(diag(outer(row(T1)[,1]+1, row(T1)[,1], "*"))),c(k-1,k))
    T1[col(T1) > row(T1)] <- 0
    T1[col(T1) == row(T1)+1] <- -(row(T1)[,1])/sqrt(diag(outer(row(T1)[,1]+1, row(T1)[,1], "*")))

    OT <- rbind(T0, T1)

    OT
}


