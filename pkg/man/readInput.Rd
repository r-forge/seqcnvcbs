\name{readInput}
\alias{readInput}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Manage reading and merging of raw datasets. Main file input
}
\description{
This is used to control the read of a meta file containing names of data files, merge, and give usable output for the main program
}
\usage{
readInput(inputFilename, formatName="Chiang", sep = "\t")
}
\arguments{
  \item{inputFilename}{
The name of file, containing relevant information of all input files
}
  \item{formatName}{
The format in which the data files are written in. We use the simple 'Chiang' as default format of input.
}
  \item{sep}{
Delimiter of the meta input file, default is tab-delimited
}
}
\details{
The meta input file should be organized in a table format with 2 columns, one of which is 'file' and the other is 'type', indicating the data file names and whether the data is from normal or tumor.
}
\value{

\item{normalSeq }{A list containing the combined normal/control reads}
\item{tumorSeq }{A list containing the combined case/tumor reads}
}
\author{
Jeremy J. Shen
}
\examples{
# This shows the format of the meta file
data(JSSim_Meta)
print(JSSim_Meta)

# This shows the recommended format, the Chiang data format
data(JSSim_NormalSim1)
print(head(JSSim_NormalSim1))
}
\seealso{
\code{\link{readListInputFile}}, \code{\link{readSeq}}
}