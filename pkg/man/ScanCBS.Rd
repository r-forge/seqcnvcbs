\name{ScanCBS}
\alias{ScanCBS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Main CBS Algorithm for Change-Point Detection
}
\description{
This is the main algorithm. It teratively scans for window of arbitrary size where the case and control read depths are different. It continues until a stopping criterion based on mBIC, maximum number of cut, and the statistic at the current segment.
}
\usage{
ScanCBS(cases, controls, statistic = "binomial", grid.size = "auto", takeN = 5, maxNCut = 30, minStat = 5, verbose = FALSE, timing = TRUE)
}
\arguments{
  \item{cases}{
A numeric vector of the case/tumor reads
}
  \item{controls}{
A numeric vector of the control/normal reads
}
  \item{statistic}{
The statistic to be used. Can be 'binomial','rabinowitz', 'normal', or 'hybrid'
}
  \item{grid.size}{
The set of grid sizes for the iterative search. An automatic default can be computed.
}
  \item{takeN}{
The number of candidate change points to be added to a temporary set at each grid size
}
  \item{maxNCut}{
The maximum number of segmentation steps to perform
}
  \item{minStat}{
The minimum statistic value required to continue the segmentation. Can be set to 0 to be ignored.
}
  \item{verbose}{
If \code{TRUE}, then will print much information on each segmentation. For diagnostics only.
}
  \item{timing}{
If \code{TRUE}, perform a timing of this algorithm, include in the output data file.
}
}
\details{
This algorithm is an use of the Circular Binary Segmentation method. It continues to segment the reads and consider the resulting child regions for further segmentation. It keeps track of the most promising cut in each children, and only the child region with the most significant segmentation is further cut, yielding more children. This is repeated until stopping criteria are met. The four types of statistics are by the use of exact binomial likelihood, score statistic, normal approximation, and an hybrid of normal approximation with binomial likelihood for small windows.
}
\value{
  \item{tauHat }{The change points called}
  \item{statHat }{A matrix containing the statistic and its segmentation for the model called}
  \item{relCN }{The relative CN computed for each segment between change points}
  \item{timingRes }{A list containing the result of the timing of this algorithm}
}
\references{
D. Rabinowitz, IMS Lecture Notes - Monograph Series, Vol. 23, 1994
}
\author{
Jeremy J. Shen
}
\seealso{
\code{\link{ScanIterateGrid}}, \code{\link{ScanBIC}}, \code{\link{relCNComp}}, \code{\link{getAutoGridSize}}
}