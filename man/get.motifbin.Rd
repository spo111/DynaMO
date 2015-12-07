\name{get.motifbin}
\alias{get.motifbin}
\title{
Divide motif sites into bins
}
\description{
Import motif sites and divide each site into consecutive bins.
}
\usage{
get.motifbin(motifseg, bin = 50)
}
\arguments{
  \item{motifseg}{
The name of a GRanges object storing motif sites. Motif sites are extended,
and usually from the output of get.motif function.
}
  \item{bin}{
The length of each bin. Default is 50 bp.
}
}
\details{
The function imports a list of motif sites stored in a GRanges object and
divides each motif site into consecutive small bins and outputs another GRanges
object.
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
motifbin1<-get.motifbin(motif1,20)
}
\keyword{ get.motifbin }
