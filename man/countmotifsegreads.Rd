\name{countmotifsegreads}
\alias{countmotifsegreads}
\title{
Calculate numbers of reads overlapping motif sites
}
\description{
Count numbers of reads overlapping motif sites.
}
\usage{
countmotifsegreads(read, segmentfile, seqL=150, core=1, format="bam",
histonemotifRFid)
}
\arguments{
  \item{read}{
A list object where each element is a GRanges object storing read coordinates
of one sample, or a data.frame object where each row is the path to a file
storing aligned reads of one sample.
}
  \item{segmentfile}{
A list object where each element is a GRanges object storing the locations of
motif sites from one motif.
}
  \item{seqL}{
The length to be extended from the 5' ends of reads.
}
  \item{core}{
The number of cores to use.
}
  \item{format}{
The format of aligned reads. The format can be a .txt file or a .bam file.
}
  \item{histonemotifRFid}{
The list object where each element contains the indexes of motif sites for
read counting. The indexes are also used to name the motif sites.
}
}
\details{
The function calculates the numbers of reads overlapping motif sites. Motif
sites from different motifs are stored in a list object. If multiple samples of
aligned reads are provided, read counts are normalized by total number of reads
across samples.
}
\value{
A list object where each element is a matrix of normalized read counts at motif
sites. Each row is the read counts at one motif site and each column is one
sample.
}
\author{
Zheng Kuang
}
\examples{
readfile<-system.file("extdata","read1.txt",package="DynaMO",mustWork=TRUE)
read1<-get.reads(readfile,150,"txt")
readfile<-system.file("extdata","motif_1.txt",package="DynaMO",mustWork=TRUE)
motif1<-get.motif(readfile,300)
motifreads<-countmotifsegreads(list(read1),list(motif1),150,1,"txt",
list(1:length(motif1)))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ countmotifsegreads }
