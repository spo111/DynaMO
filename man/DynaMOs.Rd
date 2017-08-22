\name{DynaMOs}
\alias{DynaMOs}
\title{
Predicts binding sites of transcription factors
}
\description{
DynaMOs predicts transcription factor (TF) binding sites by learning the spatial
distribution of epigenomic data with a random forest model.
}
\usage{
DynaMOs(reads, peak, motifseg, core, readsmem, markernum, time, mode, format, motifbin, readlen,
 motifname)
}
\arguments{
  \item{reads}{
A list of data.frame objects where each row is a path to a file storing aligned reads. Each data.frame
stores one chromatin mark. Or a list object with length equal to the number of histone markers.
Each element is a list object with length equal to the number of time points. Each sub-element
is a GRanges object storing locations of aligned reads.
}
  \item{peak}{
A list object with length equal to the number of histone markers. Each element
is a list object with length equal to the number of samples. Each sub-element
is a GRanges object storing locations of peaks of histone markers.
}
  \item{motifseg}{
A list object. Each element is a GRanges object storing locations of motif sites.
}
  \item{core}{
The number of cores to use.
}
  \item{readsmem}{
A logcial. If TRUE, aligned reads are imported and kept in memory. If FALSE, reads is a list of
data.frame objects.
}
  \item{markernum}{
An integer representing the number of chromatin marks.
}
  \item{time}{
An integer representing the number of time points.
}
  \item{mode}{
"s" or "l". "l" means all motif sites are examined by random forest models and
evaluated as real binding sites. "s" means only motif sites overlapping histone
marker peaks are examined and evaluated. "s" is useful when the number of motifs
and the number of motif sites are large.
}
  \item{format}{
"bam" or "txt". The format of aligned reads can be a .txt file or a .bam file.
}
  \item{motifbin}{
An integer. The length of each bin.
}
  \item{readlen}{
The length to be extended from the 5' ends of reads.
}
  \item{motifname}{
A vector of characters representing the names of each motif.
}
}
\details{
This function is part of the DynaMO pipeline. It imports location information of
peak and motif sites and aligned reads and outputs probabilities of each motif in
each sample as false discovery rates and inner products.
}
\value{
Read counts at motif sites are written in files named "motif_[motif index]_
readcount.txt". Each row is a motif site and named by the index of motif sites.
Each column is a sample.
False discovery rates are written in files named "motif_[motif index]_fdr.txt".
Rows are motif sites with the same order as readcount files and columns
represent samples.
Inner products are written in files named "motif_[motif index]_innerprod.txt".
Rows are motif sites with the same order as readcount files and columns
represent samples.
}
\author{
Zheng Kuang
}
\examples{
#Import reads, motif sites and call peaks first
library(BayesPeak)
histonelist=vector("list",length=2)
histonelist[[1]]=data.frame(c(system.file("extdata","read1.txt"
,package="DynaMO",mustWork=TRUE),system.file("extdata","read2.txt"
,package="DynaMO",mustWork=TRUE)))
histonelist[[2]]=data.frame(c(system.file("extdata","read3.txt"
,package="DynaMO",mustWork=TRUE),system.file("extdata","read4.txt"
,package="DynaMO",mustWork=TRUE)))
motif=vector(length=5)
for(i in 1:5){
	motif[i]=system.file("extdata",paste("motif_",i,".txt",sep=""
	),package="DynaMO",mustWork=TRUE)
}
motifseg=vector("list",length=motifnum)
    get.motifmulti=function(x){
        result=get.motif(motif[x],motiflen)
        return(result)
    }
motifseg=mclapply(1:motifnum,get.motifmulti,mc.cores=core)
peak=vector("list",length=2)
for(i in 1:2){
    peak[[i]]=vector("list",length=nrow(histonelist[[i]]))
    for(j in 1:nrow(histonelist[[i]])){
        tempreads=get.reads(histonelist[[i]][j,1],150,"txt")
        tempreads1=as.data.frame(tempreads)
        tempreads2=cbind(tempreads1[,1:3],tempreads1[,5])
        colnames(tempreads2)=c("chr","start","end","strand")
        temppeak=summarize.peaks(bayespeak(tempreads2))
        peak[[i]][[j]]=GRanges(seqnames=Rle(as.character(space
        (temppeak))),ranges=IRanges(start=start(temppeak),
        end=end(temppeak)))
    }
}
#DynaMOs function
DynaMOs(histonelist,peak,motifseg,2,readsmem=F,2,2,"l","txt",20,150,1:5)
}

\keyword{ DynaMOs }
