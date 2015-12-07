\name{get.reads}
\alias{get.reads}
\title{
Import aligned reads
}
\description{
Read a file containing aligned reads and convert it into a GRanges object.
}
\usage{
get.reads(file, seqL, readformat)
}
\arguments{
  \item{file}{
The name of the file where the aligned reads are stored.
}
  \item{seqL}{
The length to be extended from the 5' ends of reads.
}
  \item{readformat}{
The format of read files. The format can be a .txt file or a .bam file.
}
}
\details{
The purpose of this function is to import aligned reads and convert to GRanges
objects. Two formats can be recognized by this function, a .txt file or a .bam
file. Three columns separated by tabs are expected in a .txt file. The first
column stores chromosome name. The second column stores left end coordinates of
reads and the third column stores the strands of reads. seqL equals to the size
of DNA fragments, and usually ranges from 150 to 300.
}
\value{
A GRanges object.
}
\author{
Zheng Kuang
}

\examples{
readfile<-system.file("extdata","read1.txt",package="DynaMO",mustWork=TRUE)
read1<-get.reads(readfile,150,"txt")
}
\keyword{ get.reads }
