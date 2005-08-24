## (c) 2004-2005 Yu Chuan Tai
## University of California, Berkeley
## Created: August 2004
## Last Modified: July 4, 2005

matrix.cov <- function(x, k, trans=TRUE, c.grp=NULL, use="complete.obs")
{
   # x is already sorted by c.grp, n.grp and k.grp
   if(!is.numeric(x))  x <- as.numeric(x)
   if(missing(k)) stop("The number of time points is missing.")
   #if(k <2 & length(unique(c.grp))==1) stop("Input at least two time points.")

   #if(is.null(size)) size <- length(x)/k
   
   if(length(unique(c.grp))==1)
   {
      res <-cov(matrix(x, byrow=TRUE, ncol=k),use=use)

     if(trans)
     {
        OT1 <- ot.helmert(k)[2:k,]
        if(k==2) OT1 <- t(as.matrix(OT1))
        res <- OT1%*%res%*%t(OT1)
      }

   }

   if(length(unique(c.grp))>1 & !trans)
   {
      D <- length(unique(c.grp))
      size <- max.size <- NULL
      for(i in 1:D)
      {
         grp.indx <- c.grp==sort(unique(c.grp))[i]
         if(k==1) size[i] <- sum(apply(t(as.matrix(!apply(matrix(x[grp.indx], byrow=TRUE,ncol=k),1, is.na))),2,sum)==k)
         if(k>1) size[i] <- sum(apply((!apply(matrix(x[grp.indx], byrow=TRUE,ncol=k),1, is.na)),2,sum)==k)
         max.size[i] <- sum(grp.indx)/k
      }
      cumsize <- cumsum(max.size)
      cumsize1 <- cumsum(size)
      res1 <- as.list(NULL)
      for(i in 1:D)
      {
        if(size[i]==1) res1[[i]] <- matrix(0,k,k)
        if(i==1 & size[i]>1) res1[[i]] <- cov(matrix(x[1:(k*cumsize[1])], byrow=TRUE, ncol=k),use=use)
        if(i>1& size[i]>1) res1[[i]] <- cov(matrix(x[(k*cumsize[i-1]+1):(k*cumsize[i])], byrow=TRUE, ncol=k),use=use)
      }
      
      if(k>1) res <- matrix(apply(sapply(1:D, function(x) (size[x]-1)*res1[[x]]),1,sum),ncol=k)/(cumsize1[D]-D)
      if(k==1) res <- matrix(sum(sapply(1:D, function(x) (size[x]-1)*res1[[x]])),ncol=k)/(cumsize1[D]-D)
   }

   res

}




