\name{fitfunc}
\alias{fitfunc}
\title{
Fit a curve of a dynamic pattern
}
\description{
Take a vector of values and a vector of time and output the fitted curve, max
and max derivative time.
}
\usage{
fitfunc(expr,time)
}
\arguments{
  \item{expr}{
A vector of expression value.
}
  \item{time}{
A vector of time.
}
}
\details{
This function takes a vector of values and a vector of time and output the fitted
curve, max and max derivative time.
}
\value{
A list containing the max and max derivative time, fitted value and original value
in two data.frames with 2 columns. The first is time and the second is fitted values.
}
\author{
Zheng Kuang
}
\examples{
readfile<-system.file("extdata","motif_filter_reads_cluster.txt",package="DynaMO",mustWork=TRUE)
data<-read.table(as.character(file),sep="\t",header=T)
fitted<-fitfunc(data[1,1:4],1:4)
}
\keyword{ fitfunc }
