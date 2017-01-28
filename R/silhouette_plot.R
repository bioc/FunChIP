silhouette_plot <- function(object, p = 1, weight = NULL, 
                              alpha = 1, rescale = FALSE,
                              t.max = 0.5)
{
    if (class(object) != "GRanges")
    {
        stop('The first object is not a GRanges object.')
    }
    
    if ( (p != 0) && (p != 1) && (p != 2) )
    {
        stop ('invalid value for p, It must be 0,1,2 ')
    }
    
    if ( (alpha < 0) || (alpha > 1) )
    {
        stop ('invalid value for alpha. It must be included in [0, 1]')
    }
    
    if(is.null(weight))
    {
        stop('Provide the value for the weight used in the classification ')
    }
    
    NOregisterd <- !(is.null(object$cluster_NOshift))
    registered <- !(is.null(object$cluster_shift))
    
    maxcl <- NULL
    
    if(!(NOregisterd) && !registered)
    {
        stop('No classification provided.')
    }
    
    if(NOregisterd)
    {

        print("Silhouette for the non aligned peaks")
        dist_peaks_NO_shift <- distance_peak(object, p = p)
        cluster_NO_shift <- object$cluster_NOshift
        
        # maximum number of clusters
        maxcl <- length(cluster_NO_shift[[1]])
        if(maxcl==1)
        {
            stop('No classification provided. ')
        }
        
        dist_matrix <- (1-alpha)*dist_peaks_NO_shift$dist_matrix_d0 +
            alpha * weight * dist_peaks_NO_shift$dist_matrix_d1
        
        silh_NO_shift <- sapply(2:maxcl, 
                                function(x){compute_silhouette(x, cluster_NO_shift, dist_matrix)}, 
                                simplify = FALSE)
    }
    
    if(registered)
    {
        print("Silhouette for the aligned peaks")
        
        cluster_shift <- object$cluster_shift
        dist_shift <- object$dist_shift
        
        # maximum number of clusters
        maxcl <- length(cluster_shift[[1]])
        if(maxcl==1)
        {
            stop('No classification provided. ')
        }
        
        
        silh_shift <- sapply(2:maxcl, function(kk){
            
            
            templates <- sapply(1:kk, function(s)
            {
                dist_cluster <- sapply(dist_shift, function(y){y[s]})
                temp <- which(dist_cluster==0)
                return(temp)
            }, simplify = FALSE)
            
            
            cl_shi <- sapply(cluster_shift, function(x){x[kk]})
            define_templates <- object[templates[[kk]],]
            coefs_shift_original <- unlist(sapply(define_templates$coef_shift, function(x){x[kk]}))
            
            center_templates <- summit_peak(define_templates, 
                                            define_templates$summit_spline - coefs_shift_original)
            
            # center templates
            templates_cluster <- cluster_peak(center_templates, 
                                              parallel = FALSE, n.clust = 1, 
                                              shift.peak = TRUE, weight = weight, 
                                              alpha = alpha, p = p,
                                              t.max = t.max, plot.graph.k = FALSE, 
                                              verbose = FALSE, rescale = FALSE)
            
            coefs_shift <- unlist(templates_cluster$coef_shift)
            templates_aligned <- summit_peak(center_templates, center_templates$summit_spline - coefs_shift)
            
            # shift original peaks with the coefficient 
            # estimated by the clustering algorithm
            coef_shift_all <- sapply( object$coef_shift, function(x){x[kk]}, simplify = TRUE)
            shifted_peaks <- summit_peak(object, 
                                         object$summit_spline - coef_shift_all)
            
            # shift peaks with the coefficient 
            # related to the templates aligment
            coefs_alignment_data <- sapply(cl_shi, function(x){coefs_shift[x]})
            data_aligned <- summit_peak(shifted_peaks, 
                                        shifted_peaks$summit_spline - coefs_alignment_data)
            
            dist_peaks_shift <- distance_peak(data_aligned, p = p)
            
            dist_matrix_shift <- (1-alpha)*dist_peaks_shift$dist_matrix_d0 +
                alpha * weight * dist_peaks_shift$dist_matrix_d1
            
            shil_kk <- compute_silhouette(kk, cluster_shift, dist_matrix_shift)
            
            return(shil_kk)
        }, simplify = FALSE)
        
    }
    
    # library(RColorBrewer)
    
    if(registered && NOregisterd)
    {
        par(mfcol =c(2, maxcl-1),
            oma = c(0,0,4,0), mar=c(1.5,4,1.5,1.5))  
    }else
    {
         par(mfcol =c(1, maxcl-1),
                oma = c(0,0,4,0), mar=c(1.5,4,1.5,1.5))
    
    }
    
    for(i in 2:maxcl)
    {
        # par(oma=c(0,0,0,0))
        if(NOregisterd)
        {
            cluster_here_NO_shift <- sapply(cluster_NO_shift, function(x){x[i]})
            # plot of the silhouette index for the non aligned data
            order_data <- NULL
            shil_here <- silh_NO_shift[[i-1]]
            shil_ordered <- NULL
            for(k in 1:i)
            {
                data_here <- which(cluster_here_NO_shift == k)
                shil_ordered_here <- sort(shil_here[data_here], decreasing = TRUE)
                order_data <- c(order_data, data_here)
                shil_ordered <- c(shil_ordered, shil_ordered_here)
            }
            max_shil <- shil_ordered[c(1, (cumsum(table(cluster_here_NO_shift))[-i]+1))]
            max_shil_sort <- sort.int(max_shil, index.return = TRUE, decreasing = TRUE)
            
            new_sort <- NULL
            for(k in 1:i)
            {
                ends <- cumsum(table(cluster_here_NO_shift))
                starts <- c(1, ends[-i]+1)
                new_sort <- c(new_sort, 
                              starts[max_shil_sort$ix[k]]:ends[max_shil_sort$ix[k]])
            }
            
            scale_col <- gray.colors(i, start = 0.3, end = 0.6, gamma = 2.2, alpha = NULL)
            colors_here <- rep(scale_col, table(cluster_here_NO_shift))
            barplot(height = shil_ordered[new_sort], 
                    horiz = FALSE, 
                    col = colors_here[new_sort], border = colors_here[new_sort], 
                    ylim = c(-1, 1), cex.axis =1.5 , ylab = 'silhouette index', cex.lab = 1.5)
            title(main = paste('k =', as.character(i)) , cex.main = 2, font.main = 1) 
            abline(h = mean(shil_ordered), lwd =2 , lty = 2, col ='grey12')
            text(x = 1.1*length(shil_ordered), 
                 y = mean(shil_ordered)+ 0.08, 
                 labels = as.character(round(mean(shil_ordered),2)), cex = 1.5)
            legend('bottom', fill=scale_col, border = 0, box.lwd = 0, 
                   legend= 1:i, cex = 1.5, title = "Cluster", ncol =ceiling(i/2))
            
        }
        
        # plot aligned data
        
        if(registered)
        {
            cluster_here_shift <- sapply(cluster_shift, function(x){x[i]})
            # plot of the silhouette index for the non aligned data
            order_data <- NULL
            shil_here <- silh_shift[[i-1]]
            shil_ordered <- NULL
            for(k in 1:i)
            {
                data_here <- which(cluster_here_shift == k)
                shil_ordered_here <- sort(shil_here[data_here], decreasing = TRUE)
                order_data <- c(order_data, data_here)
                shil_ordered <- c(shil_ordered, shil_ordered_here)
            }
            max_shil <- shil_ordered[c(1, (cumsum(table(cluster_here_shift))[-i]+1))]
            max_shil_sort <- sort.int(max_shil, index.return = TRUE, decreasing = TRUE)
            
            new_sort <- NULL
            for(k in 1:i)
            {
                ends <- cumsum(table(cluster_here_shift))
                starts <- c(1, ends[-i]+1)
                new_sort <- c(new_sort, 
                              starts[max_shil_sort$ix[k]]:ends[max_shil_sort$ix[k]])
            }
            
            scale_col <- rev(rep(brewer.pal(i+3, "OrRd")[-(1:3)], length = i)[1:i])
            colors_here <- rep(scale_col, table(cluster_here_shift))
            barplot(height = shil_ordered[new_sort], 
                    horiz = FALSE, 
                    col = colors_here[new_sort], border = colors_here[new_sort], 
                    ylim = c(-1, 1), cex.axis =1.5 , ylab = 'silhouette index', cex.lab = 1.5)
            if(!NOregisterd)
            {
                title(main = paste('k =', as.character(i)) , cex.main = 2, font.main = 1) 
            }
            abline(h = mean(shil_ordered), lwd =2 , lty = 2, col ='darkred')
            text(x = 1.1*length(shil_ordered), 
                 y = mean(shil_ordered)+ 0.08, col = 'darkred',
                 labels = as.character(round(mean(shil_ordered),2)), cex = 1.5)
            legend('bottom', fill=scale_col, border = 0, box.lwd = 0, 
                   legend= 1:i, cex = 1.5, title = "Cluster", ncol = ceiling(i/2))
            
        }
        
    }
    title("Silhouette plot", outer = TRUE, font.main = 1, cex.main = 2)
    
    
    if(registered && NOregisterd)
    {
        return(list(silh_shift = silh_shift, silh_NO_shift = silh_NO_shift))
    }else
    {
        if(registered)
        {
            return(silh_shift)
        }else
        {
            return(silh_NO_shift)
        }
    }
}




###########################
### AUXILIARY FUNCTIONS ###
###########################

compute_average_distance <- function(dist, cl, cluster_i)
{
    MaxCl <- length(unique(cl))
    average_distances <- tapply(dist, cl, mean)
    ordered_aver <- c(average_distances[cluster_i], 
                      sort(average_distances[-cluster_i], decreasing = FALSE))
    return(ordered_aver)
}
compute_silhouette <- function(id_clust, cluster, dist_matrix)
{
    if(id_clust == 1)
    {
        stop('Silhouette index can not be computed with 1 cluster only.')
    }
    
    cl <- sapply(cluster, function(x){x[id_clust]})
    
    # on all the peaks compute the average distance within and
    # without the cluster
    mat_dist <- sapply(1:dim(dist_matrix)[2], 
                       function(x){compute_average_distance(
                           dist_matrix[,x], 
                           cl, cl[x])}, 
                       simplify = TRUE)
    # print(head(t(mat_dist)))
    
    s <- apply(mat_dist, 2, function(x){(x[2]-x[1])/max(x[1:2])})
    # print(head(s))
    
    return(s)
}
