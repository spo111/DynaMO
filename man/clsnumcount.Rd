\name{clsnumcount}
\alias{clsnumcount}
\title{
Calculate the number of clusters for K-means clustering
}
\description{
Take a matrix with rows representing signals at each motif site and columns
representing different samples.
}
\usage{
clsnumcount(data)
}
\arguments{
  \item{data}{
A matrix containing the signals at each motif site and sample.
}
}
\details{
This function calculates the number of clusters for K-means clustering.
}
\value{
An integer.
}
\author{
Zheng Kuang
}

\examples{
readfile<-system.file("extdata","motif_filter_reads_cluster.txt",package="DynaMO",mustWork=TRUE)
data<-read.table(as.character(file),sep="\t",header=T)
clsnum<-clsnumcount(data[,1:4])
}
\keyword{ clsnumcount }
