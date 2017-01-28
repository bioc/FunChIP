\name{bending_index}
\alias{bending_index}
\title{
The elbow rule to define the proper number of clusters.
}
\description{
Given a \linkS4class{GRanges} object with metadata columns related
to the classification performed with the \link[=cluster_peak-method]{cluster_peak} method, this function quantifies the elbow rule. See Details for a short presentation of the method and the Vignette of the package for a complete defintion of the index. }
\usage{
bending_index(object, plot.graph.k = FALSE)
}
\arguments{
  \item{object}{ \linkS4class{GRanges} object. It must contain the metadata columns associated to the classification to be analyzed. 
   Specifically it must contain the \code{cluster_NOshift} metadata (and the correspondent set of distances \code{distance_NOshift}) 
   if the user 
   wants to compute the bending index for the non aligned
   peaks and/or the \code{cluster_shift} metadata (and the correspondent set of distances \code{distance_shift}) if the user wants to compute
   the bending index for the classification with alignment.
  }
  \item{plot.graph.k}{logical. If \code{TRUE} the graph of the global distance between the
  data and corresponding center of the cluster, varying the number of clusters is plotted. Distances are normalized with the total number of peaks \eqn{n}. These are the distances uesd to compute the 
  bending index, as presented in Details.
  If \code{object} contains both the
  results with and without the classification, two lines are drown to show, beside the variation of the distance with an increase of \eqn{k}, also the decrease of the
  global distance introduced by the alignment procedure. If a single classification is stored in the \code{object}, only one line is drawn. Default is \code{FALSE}.
  }
}
\details{This function consists of the computation for each feasible value of k (from 2 to K ??? 1, with K the maximum number of clusters) of an index that quantifies the magnitude of the elbow. As higher is this index, as the correspondent value of k is meaningful. Specifically it is computed as the distance of the point in k of the global distance function (normalized with the maximum value it assumes) from the line passing by the point in k ??? 1 and in k + 1. For further details, see the Vignette.
}
\value{
The function returns 
\itemize{
    \item a data.frame (or a list with two data.frames, in case of \code{object} with  classification with and without alignment) containing the bending index for different values of the parameter \eqn{k}.
    \item if \code{plot.graph.k = TRUE} the graphical representation of the distances (normalized with the total number of peaks \eqn{n}), varying the classification type and the number of clusters.
}
}
\author{
Alice Parodi, Marco J. Morelli, Laura M. Sangalli, Piercesare Secchi, Simone Vantini
}

\examples{
# load the data
data(peaks)

# compute the bending index
index <- bending_index(peaks.data.cluster, plot.graph.k = FALSE)
# from the analysis of this results, a choice of k=3 for
# the classification with shift and k=2 for the classification 
# without shift is suggested.
}