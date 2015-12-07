\name{get.motif}
\alias{get.motif}
\title{
Read motif sites
}
\description{
Read a file containing locations of motif sites and convert it into a GRanges
object.
}
\usage{
get.motif(file, extendL)
}
\arguments{
  \item{file}{
The name of the file where the locations of motif sites are stored.
}
  \item{extendL}{
A integer. The length to be extended from the centers of motif sites.
}
}
\details{
The function reads a file containing locations of motif sites and converts it
to a GRanges object. Because real binding sites are often located near peaks of
histone modifications but not in the middle of peaks. By extending motif sites
by 200 to 1000 bp on both sides, DynaMO can capture real motif sites.
}
\value{
A GRanges object.
}
\author{
Zheng Kuang
}

\examples{
readfile<-system.file("extdata","motif_1.txt",package="DynaMO",mustWork=TRUE)
motif1<-get.motif(readfile,300)
}
\keyword{ get.motif }
