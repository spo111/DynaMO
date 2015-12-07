\name{motifbincount}
\alias{motifbincount}
\title{
Count reads overlapping consecutive bins at motif sites
}
\description{
Count numbers of reads overlapping consecutive bins at motif sites.
}
\usage{
motifbincount(motifbin, reads, core, seqL, format = "bam", name)
}
\arguments{
  \item{motifbin}{
A list object where each element is a GRanges object storing locations of
consecutive bins at motif sites of one motif.
}
  \item{reads}{
A list object where each element is a GRanges object storing read coordinates
of one sample, or a data.frame object where each row is the path to a file
storing aligned reads of one sample}
  \item{core}{
The number of cores to use.
}
  \item{seqL}{
The length to be extended from the 5' ends of reads.
}
  \item{format}{
The format of aligned reads. The format can be a .txt file or a .bam file.
}
  \item{name}{
The id of histone modification or other marker. For example, for the first
histone modification, name equals to 1. For the second histone modification,
name equals to 2.
}
}
\details{
The function calculates the numbers of reads from each sample overlapping
consecutive bins at motif sites of each motif. Two types of bin input are
expected. For the first type, each element of the input list is a GRanges
object storing bins at motif sites and the output is written in a .txt file
named motifbin_name_sampleid_motifid.txt. For the second type, each element of
the input list is still a list object containing three GRanges objects storing
positive, negative and background bins for training purpose. The outputs are
written in three .txt files names motifbin_name_sampleid_motifid_1.txt
(positive), motifbin_name_sampleid_motifid_2.txt (negative), and motifbin_name
_sampleid_motifid_3.txt (background).
}
\value{
Each output is a vector of read counts at bins of each motif from each sample
and written in a .txt file.
}
\author{
Zheng Kuang
}
\examples{
readfile<-system.file("extdata","read1.txt",package="DynaMO",mustWork=TRUE)
read1<-get.reads(readfile,150,"txt")
readfile<-system.file("extdata","motif_1.txt",package="DynaMO",mustWork=TRUE)
motif1<-get.motif(readfile,300)
motifbin1<-get.motifbin(motif1,20)
motifbincount(list(motifbin1),list(read1),1,150,"txt",1)
}
\keyword{ motifbincount }
