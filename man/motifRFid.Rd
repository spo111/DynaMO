\name{motifRFid}
\alias{motifRFid}
\title{
Find motif sites for training random forest models
}
\description{
Find ids of positive, negative and background motif sites of each motif for
training random forest models.
}
\usage{
motifRFid(motifseg, histonemotifsegreads, motifpeakoverlap, motifseg500overlapid)
}
\arguments{
  \item{motifseg}{
A list object where each element is a GRanges object storing the locations of
motif sites from one motif.
}
  \item{histonemotifsegreads}{
A list object with length equal to the number of histone markers and each
element is a list object with length equal to the number of motifs where each
sub-element is a matrix storing read counts of one histone modification at
motif sites from different motifs. Columns of matrices are read counts from
different samples and rows are different motif sites.
}
  \item{motifpeakoverlap}{
A list object with length equal to the number of histone markers and each
element contains overlapping information of motif sites and one marker. Each
element is a list object with length equal to the number of motifs and each
sub-element is a matrix containing overlapping information of motif sites from
one motif. Each column of the matrix is the numbers of peaks overlapping each
motif site.
}
  \item{motifseg500overlapid}{
A list object with length equal to the number of motifs. Each element is a
vector of numbers of peaks overlapping motif sites from this motif.
}
}
\details{
This function finds motif sites for training random forest models. It
identifies at most 250 motif sites from each motif which overlap peaks of at
least one histone marker as positive sites and 250 motif sites which don't
overlap any peak from any histone marker as negative sites. It also identifies
250 motif sites as background sites for each motif.
}
\value{
It returns a list with length equal to the number of motifs and each element is
another list with length equal to the number of samples. Each list is another
list of 3 elements. The first is the id of positive motif sites and the second
is the id of negative motif sites and the third is the id of background sites.
}
\author{
Zheng Kuang
}
\examples{
#Prepare motifseg reads
readfile<-system.file("extdata","read1.txt",package="DynaMO",mustWork=TRUE)
read1<-get.reads(readfile,150,"txt")
readfile<-system.file("extdata","motif_1.txt",package="DynaMO",mustWork=TRUE)
motif1<-get.motif(readfile,300)
motifreads<-countmotifsegreads(list(read1),list(motif1),150,1,"txt",
list(1:length(motif1)))
#Detect peaks
library(BayesPeak)
readfile<-system.file("extdata","read1.txt",package="DynaMO",mustWork=TRUE)
tempreads<-get.reads(readfile,150,"txt")
tempreads1=as.data.frame(tempreads)
tempreads2=cbind(tempreads1[,1:3],tempreads1[,5])
colnames(tempreads2)=c("chr","start","end","strand")
temppeak=summarize.peaks(bayespeak(tempreads2))
peak=GRanges(seqnames=Rle(as.character(space(temppeak))),
ranges=IRanges(start=start(temppeak),end=end(temppeak)))
#Count overlaps of motif sites and peaks
motifpeakoverlap=data.matrix(countOverlaps(motif1,peak))
#motifRFid
motifid1=motifRFid(list(motif1),list(motifreads),list(list(motifpeakoverlap)),
list(1:1000))
}
\keyword{ motifRFid }
