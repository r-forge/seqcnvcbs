\name{nhppSimConstWindow}
\alias{nhppSimConstWindow}
\title{
Simulate a Non-Homogeneous PP with constant window spike
}
\description{
Simulate non-homogeneous Poisson processes with a number of constant-widths windows of signal spike, and run CBS on them.
}
\usage{
nhppSimConstWindow(controlRates, nSpike = 25, cptLen = seq(5, 50, by = 5), nRepeat = 5, minGain = 1.5, maxGain = 5, minLoss = 0.01, maxLoss = 0.5, pGain = 0.6)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{controlRates}{
The estimated rate of nhpp for the control
}
  \item{nSpike}{
Number of signal spikes
}
  \item{cptLen}{
Constant width of the signal (CN gain or loss) region to simulate, can be a vector of different values for which to test
}
  \item{nRepeat}{
Number of total samples to generate for each choice of the width of the signal
}
  \item{minGain}{
Minimal signal gain
}
  \item{maxGain}{
Maximal signal gain
}
  \item{minLoss}{
Minimal signal loss
}
  \item{maxLoss}{
Maximal signal loss
}
  \item{pGain}{
Proportion of the signal regions that are CN gain
}
}
\details{
This function is used in conjunction with a modified, windowed rate vector to simulate non-homogeneous Poisson processes with a number of constant-widths windows of signal spike, and run CBS on them.
One should use the \code{\link{nhppRateEstimate}} function to estimate the rate of a control sample one wishes to mimic. This function randomly choose windows of a specified constant width, and spike in signals (change points) which can be either gain or loss of copy numbers.
}
\value{
Return the results of the simulation in a list
}

\author{
Jeremy J. Shen
}


\seealso{
\code{\link{nhppSimulate}}
}