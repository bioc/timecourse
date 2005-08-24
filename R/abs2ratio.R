abs2ratio <- function(x, mn, k, c.grp, reference)
{
    # x is already sorted by c.grp, n.grp and k.grp
    max.1 <- max(mn[,1])
    max.2 <- max(mn[,2])
    N <- max.1+max.2
    x1 <- x[1:(max.1*k)]
    x2 <- x[(max.1*k+1):(N*k)]
    if(reference>sort(unique(c.grp))[1])  res <- x1-x2

    if(reference==sort(unique(c.grp))[1]) res <- x2-x1

    res
}
     

#abs2ratio <-function (dummy, M, mn, k, c.grp, reference) 
#{
#    N <- max(mn[,1])+max(mn[,2])
#    x1 <- M[dummy, 1:(max(mn[,1]) * k)]
#    x2 <- M[dummy, (max(mn[,1]) * k + 1):(N * k)]
#    if (reference > sort(unique(c.grp))[1]) 
#        res <- x1 - x2
#    if (reference == sort(unique(c.grp))[1]) 
#        res <- x2 - x1
#    res
#}

