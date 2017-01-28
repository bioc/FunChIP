bending_index <- function(object, plot.graph.k = FALSE)
{
    if (class(object) != "GRanges")
    {
        stop('The first object is not a GRanges object.')
    }
    
   #  print("Quantification of the elbow rule")
    NOregisterd <- !(is.null(object$cluster_NOshift))
    registered <- !(is.null(object$cluster_shift))
    
    maxcl <- NULL
    ratio_NOreg <- NULL
    ratio_reg <- NULL
    
    if(!(NOregisterd) && !registered)
    {
        stop('No classification provided.')
    }
    
    if(NOregisterd)
    {

        # maximum number of clusters
        maxcl <- length(object$cluster_NOshift[[1]])
        
        if(maxcl<=3)
        {
            stop('To have a meaningful analysis, provide results with 
                 higher maximum number of clusters.')
        }
        
        sum_noreg <- rep(0, maxcl)
        for (i in 1:maxcl)
        {
            sum_noreg[i] <- sum(sapply(object$dist_NOshift, function(x){x[i]})) 
        }
        
        
        ratio_NOreg_df <- rep(NA, maxcl-2)
        for(i in 2:(maxcl-1))
        {
            ratio_NOreg_df[i-1] <- distance_point_line(c(i, sum_noreg[i]/max(sum_noreg)), c(i-1, sum_noreg[i-1]/max(sum_noreg)), c(i+1, sum_noreg[i+1]/max(sum_noreg)))
            
        }
       # ratio_NOreg <- diff(sum_noreg)[-length(diff(sum_noreg))]/diff(sum_noreg)[-1]
       # ratio_NOreg <- - diff(sum_noreg)[-length(sum_noreg)]/sum_noreg[-1]
      #  ratio_NOreg_df <- data.frame(ratio = t(ratio_NOreg))
       # names(ratio_NOreg_df) <- 2:(maxcl-1)
        names(ratio_NOreg_df) <- 2:(maxcl-1)
    }
    
    if(registered)
    {
        # maximum number of clusters
        maxcl <- length(object$cluster_shift[[1]])
        
        if(maxcl<=3)
        {
            stop('To have a meaningful analysis, provide results with 
                 higher maximum number of clusters.')
        }
        
        sum_reg <- rep(0, maxcl)
        for (i in 1:maxcl)
        {
            sum_reg[i] <- sum(sapply(object$dist_shift, function(x){x[i]})) 
        }
        
        #ratio_reg <- diff(sum_reg)[-length(diff(sum_reg))]/diff(sum_reg)[-1]
        #ratio_reg <- - diff(sum_reg)[-length(sum_reg)]/sum_reg[-1]
        
        ratio_reg_df <- rep(NA, maxcl-2)
        for(i in 2:(maxcl-1))
        {
            ratio_reg_df[i-1] <- distance_point_line(c(i, sum_reg[i]/max(sum_reg)), c(i-1, sum_reg[i-1]/max(sum_reg)), c(i+1, sum_reg[i+1]/max(sum_reg)))

        }
        #ratio_reg_df <- data.frame(ratio = t(ratio_reg))
        names(ratio_reg_df) <- 2:(maxcl-1)
        #names(ratio_reg_df) <- 2:(maxcl)
    }
    
    if(plot.graph.k)
    {
        if(NOregisterd && registered)
        {
            par(mar =c(5,5,4,2))
            plot(sum_reg/length(object), pch = 19, type = 'b', col = 'red3', lwd = 3, 
                 ylim=range(c(sum_reg, sum_noreg)/length(object)), 
                 ylab = 'average distance', xlab = 'number of clusters', 
                 main = 'Average distance varying k', cex = 2, 
                 cex.axis = 2, cex.lab = 2, cex.main = 2)
            lines(sum_noreg/length(object),pch = 19, 
                type = 'b', col = 'grey31', lwd = 3, cex =2)
            legend('topright', legend = c('No Shift', 'Shift'), 
                   col =c('grey31', 'red3'), lwd =3, pch = 19, 
                   cex =2, bty = 'n')
        }else
        {
            if(registered)
            {
                par(mar =c(5,5,4,2))
                plot(sum_reg/length(object), pch = 19, type = 'b', 
                     col = 'red3', lwd = 3, 
                     ylim=range(sum_reg/length(object)), 
                     ylab = 'average distance', xlab = 'number of clusters', 
                     main = 'Average distance varying k', cex = 2, 
                     cex.axis = 2, cex.lab = 2, cex.main = 2)
                legend('topright', legend = c('Shift'), 
                       col = 'red3', lwd =3, pch = 19, 
                       cex =2, bty = 'n')
            }else
            {
                par(mar =c(5,5,4,2))
                plot(sum_noreg/length(object), pch = 19, type = 'b', 
                     col = 'grey31', lwd = 3, 
                     ylim=range(sum_noreg/length(object)), 
                     ylab = 'average distance', xlab = 'number of clusters', 
                     main = 'Average distance varying k', cex = 2, 
                     cex.axis = 2, cex.lab = 2, cex.main = 2)
                legend('topright', legend = c('No Shift'), 
                       col = 'grey31', lwd =3, pch = 19, 
                       cex =2, bty = 'n')
            }
        }
        
    }
    if(registered && NOregisterd)
    {
        v <- vector('list')
        v$index_shift <- ratio_reg_df
        v$index_NOshift <- ratio_NOreg_df
       return(v) 
    }else
    {
        if(registered)
        {
            index_shift <- ratio_reg_df
            return(index_shift)
        }else
        {
            index_NOshift <- ratio_NOreg_df
            return(index_NOshift)
        }
        
    }
    
}


# AUXILIARY FUNCTION
# 
distance_point_line <- function(x, p1, p2)
{
    # check on the input points: vectors of length 2
    # (x,y)
    
    if (length(x) != 2)
    {
        stop ("x is not a point, please provide a vector 
              of length 2")
    }
    
    if (length(p1) != 2)
    {
        stop ("p1 is not a point, please provide a vector
              of length 2")
    }
    
    if (length(p2) != 2)
    {
        stop ("p2 is not a point, please provide a vector 
              of length 2")
    }
    
    num <- abs( (p2[2] - p1[2]) * x[1] -
                    (p2[1] - p1[1]) * x[2] +
                    p2[1] * p1[2] -
                    p2[2] * p1[1]
    )
    
    denom <- sqrt( (p2[2] - p1[2])^2 + 
                       (p2[1] - p1[1])^2   )
    
    return(num/denom)
}

