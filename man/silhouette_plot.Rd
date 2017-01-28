\name{silhouette_plot}
\alias{silhouette_plot}
\title{
Return the silhouette index for clustered peaks
}
\description{
It computes the silhouette index for peaks stored in a 
\linkS4class{GRanges} object and classified with the \link[=cluster_peak-method]{cluster_peak} method. If the two classifications with a and without 
alignment are provided, this method computes the index for 
both these classifications.
}
\usage{
silhouette_plot(object, p = 1, 
        weight = NULL, alpha = 1, 
        rescale = FALSE, t.max = 0.5)
}
\arguments{
   \item{object}{ \linkS4class{GRanges} object. It must contain the metadata columns associated to the classification to be analyzed. 
   Specifically it must contain the \code{cluster_NOshift} metadata 
   if the user 
   wants to compute the silhouette index for the non aligned
   peaks and/or the \code{cluster_shift} metadata if the user wants to compute
   the index for the classification with alignment.
  }
  
  \item{p}{integer value in \{0, 1 , 2\}. Order of the \eqn{L^p} distance
  used. In particular \code{p = 0} stands for the \eqn{L^{\infty}} distance,
  \code{p = 1} for \eqn{L^1} and \code{p = 2} for \eqn{L^2}. Default is 1.}
  
    \item{weight}{real. Weight \eqn{w} of the distance function (see Details for the definition of the distance function), needed to
  make the distance between splines and  derivatives comparable. It has no Default since it must be the same weight used to define the distance for the classification. 
}

 \item{alpha}{real value between 0 and 1. Value of the convex weight \eqn{\alpha} of the distance to
  balance the distance between data and derivatives.
  See details for the definition. Default is 1.
  }

 \item{t.max}{real value. It tunes the maximum shift allowed. In
  particular the maximum shift at each iteration is computed as
  \deqn{
    \textrm{max\_shift} = \textrm{\code{t.max}} * \textrm{range}({\textrm{\code{object}}})
  }{max_shift = t.max * range(object)}
  and the optimum registration coefficient will be identified between - max_shift and
  + max_shift. range({\code{object}}) is the maximum amplitude of the peaks. Default is 0.5.
  }
  
  \item{rescale}{logical. If \code{TRUE} clustering is performed on \emph{scaled peaks}. For the 
  definition of \emph{scaled peaks} see \link[=smooth_peak-method]{smooth_peak}.}
}
\details{
See [Rousseeuw, 1987] for the detailed definition of the index. Specifically, for the peak \eqn{i} it is computed as
\deqn{
s(i) = \frac{a(i)-b(i)}{\max(a(i), b(i))}
}
with \eqn{a(i)} the average dissimilarity of peak \eqn{i} with all other data within the same cluster and \eqn{b(i)} the lowest average dissimilarity of  \eqn{i} to any other cluster, of which \eqn{i} is not a member.
}
\value{
The function returns 
\itemize{
    \item the list of the silhouette indeces for the two classifications (if provided in the \linkS4class{GRanges} object) and for all the choices of the number of clusters 
    \item the graphical representation of the silhouette index, varying the number of clusters and the classification. The average silhouette index is also presented.

}
}
\references{
Peter J. Rousseeuw (1987). Silhouettes: a Graphical Aid to the Interpretation and Validation of Cluster Analysis. Computational and Applied Mathematics. 20: 53???65.
}
\author{
Alice Parodi, Marco J. Morelli, Laura M. Sangalli, Piercesare Secchi, Simone Vantini
}
\examples{
# load the data
data(peaks)

# computes the silhouette index and 
# shows the graph
sil <- silhouette_plot(peaks.data.cluster, p=2, weight = 1, alpha = 1,
                         rescale = FALSE, t.max = 2)
}