## (c) 2004-2005 Yu Chuan Tai
## University of California, Berkeley
## Plots the temporal profile for a single gene
## Created: August, 2004
## Last Modified: July 23, 2005

plotProfile <- 
function (object, stats = c("HotellingT2", "MB"), ranking = 1, 
    gid = NULL, gnames = NULL, desc = NULL, type = c("p", "l", 
        "b"), col = 2:100, lty = 1:100, pch = 1:100, 
    lwd = 2, xlab = "Time", ylab = "Expression", legloc = NULL, 
    xlim = NULL, ylim = NULL, cex.main = 1, ...) 
{
    if (!is.null(object$geneNames)) 
        gnames <- object$geneNames
    if (!is.null(object$descriptions)) 
        desc <- object$descriptions
    indx <- order(object$con.group, object$rep.group, object$time.group)
    M <- object$M[, indx]
    G <- nrow(M)
    if (!is.character(desc)) 
        desc <- as.character(desc)
    stats <- match.arg(stats, c("HotellingT2", "MB"))
    type <- match.arg(type, c("p", "l", "b"))
    if (!is.null(gid)) 
        ranking <- NULL
    if (stats == "HotellingT2" & !is.null(ranking)) {
        pos <- object$pos.HotellingT2[ranking]
        val <- object$HotellingT2[pos]
        gid <- gnames[pos]
    }
    if (stats == "MB" & !is.null(ranking)) {
        pos <- object$pos.MB[ranking]
        val <- object$MB[pos]
        gid <- gnames[pos]
    }
    if (stats == "HotellingT2" & is.null(ranking)) {
        pos <- c(1:G)[gnames == gid]
        if (length(pos) > 1) 
        {
            pos <- pos[1]
            warning(length(pos), " genes have the same gene ID...only the first gene is plotted!") 
        }
        val <- object$HotellingT2[pos]
        ranking <- G - rank(object$HotellingT2)[pos] + 1
    }
    if (stats == "MB" & is.null(ranking)) {
        pos <- c(1:G)[gnames == gid]
        if (length(pos) > 1) 
        {
            pos <- pos[1]
            warning(length(pos), " genes have the same gene ID...only the first gene is plotted!") 
        }
        val <- object$MB[pos]
        ranking <- G - rank(object$MB)[pos] + 1
    }
    x <- M[pos, ]
    desc <- desc[pos]
    size <- as.matrix(object$size)[pos, ]
    D <- length(size)
    k <- length(unique(object$time.group))
    y <- matrix(x, byrow = TRUE, ncol = k)
    yulim <- max(y, na.rm = TRUE) + 0.4 * (max(y, na.rm = TRUE) - 
        min(y, na.rm = TRUE))
    yllim <- min(y, na.rm = TRUE) - 0.4 * (max(y, na.rm = TRUE) - 
        min(y, na.rm = TRUE))
    xulim <- max(unique(object$time.group), na.rm = TRUE) + 1
    xllim <- min(unique(object$time.group), na.rm = TRUE) - 1
    if (is.null(ylim)) 
        ylim <- c(yllim, yulim)
    if (is.null(xlim)) 
        xlim <- c(xllim, xulim)
    if (!is.null(desc)) {
        plot(sort(unique(object$time.group)), y[1, ], type = "n", 
            main = paste(gid, "", stats, "=", round(val, 1), 
                "", "rank=", ranking), xlim = xlim, ylim = ylim, 
            xlab = xlab, ylab = ylab, cex.main = cex.main)
        mtext(substr(desc, 1, 10), 4)
    }
    if (is.null(desc)) 
        plot(sort(unique(object$time.group)), y[1, ], type = "n", 
            main = paste(gid, "", stats, "=", round(val, 1), 
                "", "rank=", ranking), xlim = xlim, ylim = ylim, 
            xlab = xlab, ylab = ylab, cex.main = cex.main)
    rep.pch <- rep.lty <- con <- matrix(rep(0, max(as.matrix(object$size)) * 
    ncol(as.matrix(object$size))), ncol = ncol(as.matrix(object$size)))
    
    con2 <- NULL
    if (D !=2) {
        rep.indx <- unique(cbind(object$con.group, object$rep.group))
        #tmp <- (rep.indx[, 1] == sort(unique(object$con.group))[1])
        tmp <- sapply(1:D, function(z) rep.indx[,1]==sort(unique(object$con.group))[z])
        for(z in 1:D)
        {
           con <- con1 <- sort(rep.indx[tmp[,z], 2])
           con2[z] <- length(con)
           if(z==1)
           {
              rep.pch[1:con2[z], z] <- pch[1:cumsum(con2)[z]]
              rep.lty[1:con2[z], z] <- lty[1:cumsum(con2)[z]]

           }
           if(z >1 )
           {
               rep.pch[1:con2[z], z] <- pch[(cumsum(con2)[z-1]+1):cumsum(con2)[z]]
               rep.lty[1:con2[z], z] <- lty[(cumsum(con2)[z-1]+1):cumsum(con2)[z]]
            }
        }
    }

    if (D == 2) {
        rep.indx <- unique(cbind(object$con.group, object$rep.group))
        tmp <- (rep.indx[, 1] == sort(unique(object$con.group))[1])
        con1 <- sort(rep.indx[tmp, 2])
        tmp <- (rep.indx[, 1] == sort(unique(object$con.group))[2])
        con2 <- sort(rep.indx[tmp, 2])
        if (length(con1) < length(con2)) 
            con1[(length(con1) + 1):length(con2)] <- NA
        if (length(con2) < length(con1)) 
            con2[(length(con2) + 1):length(con1)] <- NA
        con <- data.frame(con1, con2)
        rep.pch[1:nrow(con[!is.na(con[, 1]), ]), 1] <- pch[1:nrow(con[!is.na(con[, 
            1]), ])]
        rep.lty[1:nrow(con[!is.na(con[, 1]), ]), 1] <- lty[1:nrow(con[!is.na(con[, 
            1]), ])]
        both <- intersect(con[, 1], con[, 2])
        if (length(both) > 0) {
            for (i in 1:length(both)) {
                for (j in 1:nrow(con)) {
                  if (con[j, 2] == both[i] & !is.na(con[j, 2])) 
                    rep.pch[j, 2] <- rep.pch[as.numeric(con[, 
                      1]) == as.numeric(both[i]), 1]
                  if (con[j, 2] == both[i] & !is.na(con[j, 2])) 
                    rep.lty[j, 2] <- rep.lty[as.numeric(con[, 
                      1]) == as.numeric(both[i]), 1]
                }
            }
        }
        for (j in 1:nrow(con)) {
            if (rep.pch[j, 2] == 0) 
                rep.pch[j, 2] <- pch[nrow(con[!is.na(con[, 1]), 
                  ]) + j]
            if (rep.lty[j, 2] == 0) 
                rep.lty[j, 2] <- lty[nrow(con[!is.na(con[, 1]), 
                  ]) + j]
        }
    }
    tmp <- NULL
    con.indx <- cumsum(apply(as.matrix(object$size),2,max))
    for (j in 1:D) {
        if (j == 1) 
            mydata <- y[1:con.indx[j], ]
        if (j > 1) 
            mydata <- y[(con.indx[j - 1] + 1):con.indx[j], ]
        #if (is.null(object$weights)) 
        #    mean.vec <- apply(mydata, 2, mean, na.rm = TRUE)
        #if (!is.null(object$weights)) {
        #    w1 <- matrix(object$weights[pos, 1:(k * size[j])], 
        #        byrow = FALSE, ncol = size[j])
        #    totalw1 <- apply(w1, 1, sum)
        #    object$weights[pos, 1:(k * size[j])] <- as.numeric(w1/totalw1)
        #    mean.vec <- apply(matrix(object$weights[pos, 1:(k * 
        #        size[j])] * as.numeric(t(mydata)), byrow = FALSE, 
        #        ncol = size[j]), 1, sum, na.rm = TRUE)
        #}
        for (i in 1:nrow(mydata)) {
            if (length(size) == 1) {
                lines(sort(unique(object$time.group)), mydata[i, 
                  ], type = type, col = col[i], lty = rep.lty[i, 
                  j], pch = rep.pch[i, j], lwd = lwd)
                tmp <- cbind(tmp, sort(unique(object$time.group)), 
                  mydata[i, ])
            }
            if (D > 1) {
                lines(sort(unique(object$time.group)), mydata[i, 
                  ], type = type, col = col[j], lty = rep.lty[i, 
                  j], pch = rep.pch[i, j], lwd = lwd)
                tmp <- cbind(tmp, sort(unique(object$time.group)), 
                  mydata[i, ])
            }
        }
        #if (avg) {
        #    if (length(size) == 1) 
        #        lines(sort(unique(object$time.group)), mean.vec, 
        #          type = "l", lwd = lwd, col = 1)
        #    if (D > 1) 
        #        lines(sort(unique(object$time.group)), mean.vec, 
        #          type = "l", lwd = lwd, col = col[j])
        #}
    }
    corner <- matrix(c(xlim[1], ylim[1], xlim[1], ylim[2], xlim[2], 
        ylim[1], xlim[2], ylim[2]), byrow = TRUE, ncol = 2)
    tmp <- matrix(as.numeric(t(tmp)), byrow = TRUE, ncol = 2)
    loc <- NULL
    for (m in 1:nrow(corner)) {
        tt <- rbind(corner[m, ], tmp)
        corner.dist <- dist(tt)
        loc[m] <- min(as.matrix(corner.dist)[-1, 1])
    }
    if (is.null(legloc)) {
        legloc <- corner[order(loc, decreasing = TRUE)[1], ]
        if (order(loc, decreasing = TRUE)[1] == 1) 
            legloc[2] <- legloc[2] + loc[1]
        if (order(loc, decreasing = TRUE)[1] == 2) 
            legloc[2] <- legloc[2] - loc[2]
        if (order(loc, decreasing = TRUE)[1] == 3) {
            legloc[1] <- legloc[1] - loc[3]
            legloc[2] <- legloc[2] + loc[3]
        }
        if (order(loc, decreasing = TRUE)[1] == 4) 
            legloc[1] <- legloc[1] - loc[4]
    }
    if (type == "b") {
        rep.pch <- as.vector(rep.pch)[as.vector(rep.pch) > 
            0]
        rep.lty <- as.vector(rep.lty)[as.vector(rep.lty) > 
            0]
         
        if (D == 1) {
            na.indx <- apply(y,1,is.na)[1,]
            rep.lty <- rep.lty[1:max(object$size)][!na.indx]
            rep.pch <- rep.pch[1:max(object$size)][!na.indx]
            col <- col[1:max(object$size)][!na.indx]
            legend(legloc[1], legloc[2], legend = sort(unique(paste(object$con.group, 
                object$rep.group)))[!na.indx], lty = rep.lty, lwd = lwd, 
                col = col, pch = rep.pch, text.col = col, 
                cex = 0.75)
        }
        if (D > 1) {
            na.indx <- apply(y,1,is.na)[1,]
            tmp <- sum(apply(object$size,2,max))
            rep.lty <- rep.lty[1:tmp][!na.indx]
            rep.pch <- rep.pch[1:tmp][!na.indx]
            col <- col[rep(1:length(unique(object$con.group)), apply(object$size,2,max))][!na.indx]
            legend(legloc[1], legloc[2], legend = sort(unique(paste(object$con.group, 
                object$rep.group)))[!na.indx], lty = rep.lty, lwd = lwd, 
                col = col, pch = rep.pch, text.col = col, cex = 0.75)
        }
    }
    if (type == "p") {
        rep.pch <- as.vector(rep.pch)[as.vector(rep.pch) > 
            0]
        rep.lty <- as.vector(rep.lty)[as.vector(rep.lty) > 
            0]
        if (D == 1){
            na.indx <- apply(y,1,is.na)[1,]
            rep.pch <- rep.pch[1:max(object$size)][!na.indx]
            col <- col[1:max(object$size)][!na.indx]
            legend(legloc[1], legloc[2], legend = sort(unique(paste(object$con.group, 
                object$rep.group)))[!na.indx], col = col, pch = rep.pch, 
                text.col = col, cex = 0.75)
        }
        if (D > 1){
            na.indx <- apply(y,1,is.na)[1,]
            tmp <- sum(apply(object$size,2,max))
            rep.pch <- rep.pch[1:tmp][!na.indx]
            col <- col[rep(1:length(unique(object$con.group)), apply(object$size,2,max))][!na.indx]
            legend(legloc[1], legloc[2], legend = sort(unique(paste(object$con.group, 
                object$rep.group)))[!na.indx], col = col, pch = rep.pch, 
                text.col = col, cex = 0.75)
        }
    }
    if (type == "l") {
        rep.pch <- as.vector(rep.pch)[as.vector(rep.pch) > 
            0]
        rep.lty <- as.vector(rep.lty)[as.vector(rep.lty) > 
            0]
        if (D == 1) {
            na.indx <- apply(y,1,is.na)[1,]
            rep.lty <- rep.lty[1:max(object$size)][!na.indx]
            col <- col[1:max(object$size)][!na.indx]
            legend(legloc[1], legloc[2], legend = sort(unique(paste(object$con.group, 
                object$rep.group)))[!na.indx], lty = rep.lty, lwd = lwd, 
                col = col, text.col = col, cex = 0.75)
        }
        if (D > 1) {
            na.indx <- apply(y,1,is.na)[1,]
            tmp <- sum(apply(object$size,2,max))
            rep.lty <- rep.lty[1:tmp][!na.indx]
            col <- col[rep(1:length(unique(object$con.group)), apply(object$size,2,max))][!na.indx]
            legend(legloc[1], legloc[2], legend = sort(unique(paste(object$con.group, 
                object$rep.group)))[!na.indx], lty = rep.lty, lwd = lwd, 
                col = col, text.col = col, cex = 0.75)
        }
    }
}

