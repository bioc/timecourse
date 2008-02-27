## (c) 2004-2005 Yu Chuan Tai
## University of California, Berkeley
## Created: August 2004
## Last Modified: July 16, 2005

    
mb.paired <- function(object, k, mn, nu=NULL, Lambda=NULL, eta=NULL,
c.grp=NULL, k.grp=NULL, n.grp=NULL, one.D=FALSE, reference=NULL, r=FALSE, vec=FALSE, d=NULL, prop=0.02, 
T2.only=TRUE)
{

    res <- as.list(NULL)
    M <- NULL
    geneNames <- desc <- NULL
     if (is(object, "MAList") || is(object, "list"))
          M <- object$M

    if (is(object, "marrayNorm"))
    {
          M <- object@maM
          geneNames <- rownames(M)             
          desc <- object@maGnames@maLabels
    }

    if (is(object, "ExpressionSet"))
    {
          M <- exprs(object)
          geneNames <- rownames(M)       
    }

     if (is.null(M))
        M <- as.matrix(object)

    G <- nrow(M)
    if(missing(k)) stop("The number of time points is missing.")
    if(missing(mn)) stop("The sample sizes are missing.")
    if(one.D) 
    {
       max.mn <- ncol(M)/k
       if(!is.null(mn)) mn <- as.matrix(mn)
       if(!is.null(mn) & nrow(mn) !=G) stop("The sample sizes are incorrect.")

       if(is.null(c.grp)) 
       {
            c.grp <- rep(1, ncol(M))
            cat("Condition group assignments are set to default.","\n")
       }
       if(length(unique(c.grp)) != 1) stop("The biological condition group assignments are incorrect!")
       if(is.null(k.grp)) 
       {
            k.grp <- rep(1:k, max.mn)
            cat("Time group assignments are set to default.","\n")
       }
       if(length(unique(k.grp)) != k) 
            stop("Error in the number of time points or the time group assignments!")
       if(is.null(n.grp)) 
       {    
            n.grp <- rep(c(1:max.mn), rep(k, max.mn))
            cat("Replicate group assignments are set to default.","\n")
       }
       if(length(unique(n.grp)) != max.mn) 
             stop("Error in the sample size or the replicate group assignments.")
     }

    if(!one.D) 
    {
       max.mn <- ncol(M)/k/2
       #if((is.null(mn) || sum(is.na(mn))>0) & is.null(c.grp)) stop("The sample sizes are missing!")
       
       if(!is.null(mn) & !is.null(c.grp)) 
       {
          if(length(unique(c.grp)) != 2) stop("Error in biological condition group assignments!") 
          if(ncol(mn) != length(unique(c.grp))) stop("Error in sample sizes!")
          if(sum(mn[,1]!=mn[,2])>0) stop("Some genes have unequal sample sizes.")
          if(((max.mn*k) != sum(c.grp==unique(c.grp)[1])) || 
                         ((max.mn*k) != sum(c.grp==unique(c.grp)[2])))
                   stop("Error in sample sizes or biological condition group assignments!")
       } 
     
       if(is.null(c.grp)) 
       {
            c.grp <- rep(c(1,2), each=max.mn*k)
            cat("Condition group assignments are set to default.","\n")
       }
       if(length(unique(c.grp)) != 2) stop("The biological condition group assignments are incorrect!")
       if(is.null(k.grp)) 
       {
             k.grp <- rep(1:k, ncol(M)/k)
             cat("Time group assignments are set to default.","\n")
       }
       if(length(unique(k.grp)) != k) 
             stop("Error in the number of time points or the time group assignments!")
       if(is.null(n.grp)) 
       {
             n.grp <- rep(c(1:max.mn,1:max.mn), rep(k, ncol(M)/k))
             cat("Replicate group assignments are set to default.","\n")
       }
       if(length(unique(n.grp)) != max.mn) 
             stop("Error in the sample sizes or the replicate group assignments!")
       if(is.null(reference)) reference <- c.grp[1]
     }
   
   mydata <- M
   indx <- order(c.grp, n.grp, k.grp)
   M <- M[,indx]

   mis <- apply(!apply(M, 1, is.na), 2, sum)
   mis <- sum((mis/k-floor(mis/k)) !=0)
   if(mis>0) stop(mis, " genes may have missing values in some of the replicates.")


   if(!one.D)  M <- t(apply(M, 1, abs2ratio, mn, k, c.grp, reference))
   
     S <- apply(M, 1, matrix.cov, k, trans=FALSE, c.grp=rep(1,ncol(M)))
     diagS <- apply(S, 2, function(x) diag(matrix(x,ncol=k))) 
     S.avg <- matrix(apply(S, 1, mean, na.rm=TRUE),ncol=k)
  

   if(is.null(nu)||is.null(Lambda))
   {

   nu.lim <- k+6
   
   if(!is.null(nu)) 
   {
     nu0 <- nu
     nu <- max(nu0, nu.lim)
     if(is.infinite(nu) & is.null(Lambda) )  {Lambda <- S.avg} 
     if(is.finite(nu)& is.null(Lambda))    {Lambda <- (nu-k-1)*S.avg/nu}
     nu<- nu0
   }

   if(is.null(nu))
   {
     nu0 <- mean(sapply(1:k, function(x) squeezeVar(diagS[x,], mn[,1]-1)$df.prior))
     nu <- max(nu0, nu.lim)
     if(is.infinite(nu)& is.null(Lambda))  {Lambda <- S.avg} 
     if(is.finite(nu)& is.null(Lambda))    {Lambda <- (nu-k-1)*S.avg/nu}
     nu <- nu0
   }
   }
   max.n <- max(mn[,1])
   xbar <- apply(M, 1, function(x) apply(matrix(x,  byrow=FALSE, ncol=max.n),1,mean,na.rm=TRUE))
   if(r)
   {
        e <- sapply(1:G, function(x) matrix(M[x,], byrow=FALSE, ncol=max.n)-xbar[,x])
        tune <- sapply(1:G, function(x) d * median(abs(e[, x]),na.rm = TRUE)/0.6745)
        indx <- sapply(1:G, function(x) abs(e[, x]) > tune[x])
        na.indx <- sapply(1:G, function(x) c(1:(k * max.n))[!is.na(indx[,x])])
        wt <- matrix(rep(1, G * max.n * k), ncol = max.n * k)
        wt <- t(sapply(1:G, function(x) {wt[x,][-na.indx[,x]] <- 0
              wt[x,]}))
        wt <- t(sapply(1:G, function(x) {wt[x, ][na.indx[,x]][indx[,x][na.indx[,x]]] <-
        tune[x]/abs(e[,x][na.indx[,x]][indx[,x][na.indx[,x]]])
        wt[x,]}))
        totalw <- sapply(1:G, function(x) apply(matrix(wt[x,], byrow=FALSE, ncol=max.n), 1, sum,
        na.rm=TRUE))
        w <- sapply(1:G, function(x) matrix(wt[x,], byrow=FALSE, ncol=max.n)/totalw[,x])
        xbar <- sapply(1:G, function(x)
        apply(matrix(w[, x] * M[x, ], byrow = FALSE, ncol = max.n),1, sum,na.rm=TRUE))



    }


   tol <- .Machine$double.eps
   if(is.finite(nu) & nu > tol) 
   {
     modS <- sapply(1:G, function(x) ((mn[x,1]-1)*matrix(S[,x],ncol=k)+nu*Lambda)/(mn[x,1]-1+nu))
     if(is.null(eta) & (T2.only==FALSE))
     {
     sqrt.modS <- apply(modS, 2, function(x) 
     svd(matrix(x,ncol=k))$u%*%diag(sqrt(svd(matrix(x,ncol=k))$d))%*%t(svd(matrix(x,ncol=k))$v)) 
     modt <- sapply(1:G, function(x) mn[x,1]^(1/2)*solve(matrix(sqrt.modS[,x],ncol=k))%*%xbar[,x])
     HotellingT2 <- apply(modt, 2, function(x) t(x)%*%x)
     }
     if(T2.only)
     HotellingT2 <- sapply(1:G, function(x) mn[x,1]*t(xbar[,x])%*%solve(matrix(modS[,x],ncol=k))%*%xbar[,x])

     if(!is.null(eta) & (T2.only==FALSE))
         HotellingT2 <- sapply(1:G, function(x) mn[x,1]*t(xbar[,x])%*%solve(matrix(modS[,x],ncol=k))%*%xbar[,x])
   }
   if(is.infinite(nu))
   {
     modS <- sapply(1:G, function(x) Lambda)
     if(is.null(eta) & (T2.only==FALSE))
     {
     sqrt.modS <- apply(modS, 2, function(x)
     svd(matrix(x,ncol=k))$u%*%diag(sqrt(svd(matrix(x,ncol=k))$d))%*%t(svd(matrix(x,ncol=k))$v))
     modt <- sapply(1:G, function(x) mn[x,1]^(1/2)*solve(matrix(sqrt.modS[,x],ncol=k))%*%xbar[,x])
     HotellingT2 <- apply(modt, 2, function(x) t(x)%*%x)
     }
     if(T2.only) 
     HotellingT2 <- sapply(1:G, function(x) 
        mn[x,1]*t(xbar[,x])%*%solve(matrix(modS[,x],ncol=k))%*%xbar[,x])

     if(!is.null(eta) & (T2.only==FALSE))
         HotellingT2 <- sapply(1:G, function(x) mn[x,1]*t(xbar[,x])%*%solve(matrix(modS[,x],ncol=k))%*%xbar[,x])
   }
   if(nu < tol & nu >=0)           
   {
     modS <- S
     if(is.null(eta) & (T2.only==FALSE))
     {                                                     
     sqrt.modS <- apply(modS, 2, function(x)
     svd(matrix(x,ncol=k))$u%*%diag(sqrt(svd(matrix(x,ncol=k))$d))%*%t(svd(matrix(x,ncol=k))$v))
     modt <- sapply(1:G, function(x) mn[x,1]^(1/2)*ginv(matrix(sqrt.modS[,x],ncol=k))%*%xbar[,x])
     HotellingT2 <- apply(modt, 2, function(x) t(x)%*%x)
     } 
     if(T2.only)
     HotellingT2 <- sapply(1:G, function(x)
     mn[x,1]*t(xbar[,x])%*%ginv(matrix(modS[,x],ncol=k))%*%xbar[,x])

     if(!is.null(eta) & (T2.only==FALSE))
         HotellingT2 <- sapply(1:G, function(x) mn[x,1]*t(xbar[,x])%*%ginv(matrix(modS[,x],ncol=k))%*%xbar[,x])

   }

   # HotellingT2 does not require eta
   if(is.null(eta) & (T2.only==FALSE)) eta <- mean(sapply(1:k, function(x) tmixture.matrix(modt[x,], 
   sqrt(mn[,1])^(-1), mn[,1]+nu-1, prop, c(0.1, 4)^2/Lambda[x,x]))^(-1)) 
   res$M <- mydata
   res$prop <- prop
   res$nu <- nu
   res$Lambda <- Lambda
   res$eta <- eta

 
  if(nu < 0) stop("The estimation of prior degrees of freedom <0 !")
   
   if(is.finite(nu) & nu > tol)
   {
      MB <- NULL
      if(T2.only==FALSE)
      {
      MB1 <- log(prop,10)-log(1-prop,10)
      MB2 <- (k/2)*log((eta/(mn[,1]+eta)),10)
      MB3 <- ((mn[,1]+nu)/2)*log(((mn[,1]-1+nu+HotellingT2)/(mn[,1]-1+nu+(eta/(mn[,1]+eta))*HotellingT2)), 10)
      MB <- MB1+MB2+MB3
      }
      unique.mn <- unique(mn[,1])
      res$percent <- matrix(c(unique.mn, round(100*nu/(nu+unique.mn-1))), ncol=2)
      res$size <- mn
      res$con.group <- c.grp
      res$rep.group <- n.grp
      res$time.group <- k.grp
      if(r) res$weights <- wt
      if(vec) res$modt <- modt
      res$HotellingT2 <- HotellingT2 
      res$MB <- MB
   }

   if(is.infinite(nu))
   {
      MB <- NULL
      if(T2.only==FALSE)
      {
      MB1 <- log(prop,10)-log(1-prop,10)
      MB2 <- (k/2)*log((eta/(mn[,1]+eta)),10)
      MB3 <- 0.5*(mn[,1]/(mn[,1]+eta))*HotellingT2-log(10)
      MB <- MB1+MB2+MB3
      }
      unique.mn <- unique(mn[,1])
      res$percent <- matrix(c(unique.mn, 100), ncol=2)
      res$size <- mn
      res$con.group <- c.grp
      res$rep.group <- n.grp
      res$time.group <- k.grp
      if(r) res$weights <- wt
      if(vec) res$modt <- modt  
      res$HotellingT2 <- HotellingT2 
      res$MB <- MB
   }

   if(nu < tol & nu >=0)
   {
       MB <- NULL
       if(T2.only==FALSE)
       {
       MB1 <- log(prop,10)-log(1-prop,10)
       MB2 <- (k/2)*log((eta/(mn[,1]+eta)),10)
       MB3 <- (mn[,1]/2)*log(((mn[,1]-1+HotellingT2)/(mn[,1]-1+(eta/(mn[,1]+eta))*HotellingT2)), 10)
       MB <- MB1+MB2+MB3
       }
      unique.mn <- unique(mn[,1])
      res$percent <- matrix(c(unique.mn, 0), ncol=2)
      res$size <- mn
      res$con.group <- c.grp
      res$rep.group <- n.grp
      res$time.group <- k.grp
      if(r) res$weights <- wt
      if(vec) res$modt <- modt 
      res$HotellingT2 <- HotellingT2 
      res$MB <- MB
   }

   res$pos.HotellingT2 <- order(HotellingT2, decreasing = TRUE)
   if(!is.null(MB)) res$pos.MB <- order(MB, decreasing = TRUE)
   res$reference <- reference
   res$geneNames <- geneNames
   res$descriptions <- desc
   new("MArrayTC", res)

}
