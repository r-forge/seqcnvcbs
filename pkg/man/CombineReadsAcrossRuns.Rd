\name{CombineReadsAcrossRuns}
\alias{CombineReadsAcrossRuns}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Combine multiple read lists
}
\description{
Combines multiple lists in the same format of the same sample into one list of the said format.
}
\usage{
CombineReadsAcrossRuns(seqs)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{seqs}{
	A list of lists, each containing equal number of numeric vectors that can be concatenated together. Both number of lists and number of variables can be arbitrary.
}
}
\value{
Returns a list of the same format as the input lists
}
\author{
Jeremy J. Shen
}
\seealso{
\code{\link{ScanCBS}}
}