univ.func <-
function (dummy, M, k, n, indx = 1) 
{
    #if (!is.numeric(M[dummy,])) 
        x <- as.numeric(M[dummy,])
    if (missing(n)) 
        n <- apply(M, 1, function(x) sum(apply((!apply(matrix(x, byrow=TRUE, ncol=k),1, 
        is.na)),2,sum)>=1))

    x <- matrix(x, byrow = FALSE, ncol = max(n))
    OT <- ot.helmert(k)
    res <- OT[indx, ] %*% x
    res
}

