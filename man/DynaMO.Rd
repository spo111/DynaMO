\name{DynaMO}
\alias{DynaMO}
\title{
Predicts spatiotemporal binding activities of transcription factors
}
\description{
DynaMO predicts transcription factor (TF) binding sites by learning the spatial
distribution of epigenomic data with a random forest model. In a time-course
experiment, DynaMO also predicts the temporal binding patterns of TFs at
predicted binding sites. DynaMO identifies important TFs in a dynamic process.
}
\usage{
DynaMO(readlist, peak, motif, mode="l", core=1, readsmem=TRUE, readlen=150,
readformat="bam", motiflen=300, motifbin=20, cluster = 0, fdrcut=0.01,batch=T)
}
\arguments{
  \item{readlist}{
A data.frame object where each row is a path to a file storing aligned reads.
}
  \item{peak}{
A list object with length equal to the number of histone markers. Each element
is a list object with length equal to the number of samples. Each sub-element
is a GRanges object storing locations of peaks of histone markers.
}
  \item{motif}{
A vector object with length equal to the number of motifs. Each element is a
path to a file storing locations of motif sites.
}
  \item{mode}{
"s" or "l". "l" means all motif sites are examined by random forest models and
evaluated as real binding sites. "s" means only motif sites overlapping histone
marker peaks are examined and evaluated. "s" is useful when the number of motifs
and the number of motif sites are large.
}
  \item{core}{
The number of cores to use.
}
  \item{readsmem}{
A logcial. If TRUE, aligned reads are imported once and kept in memory. It
requires large memory but saves time.
}
  \item{readlen}{
The length to be extended from the 5' ends of reads.
}
  \item{readformat}{
"bam" or "txt". The format of aligned reads can be a .txt file or a .bam file.
}
  \item{motiflen}{
An integer. The length to be extended from the centers of motif sites.
}
  \item{motifbin}{
An integer. The length of each bin.
}
  \item{cluster}{
An integer. The number of clusters for K-means clustering. 0 means automatic
determination of number of clusters.
}
  \item{fdrcut}{
A numeric value. Motif sites with adjusted p values lower than fdrcut are used
for clustering.
}
  \item{batch}{
A logcial. If TRUE, motifs are divided into batches for processing. It is recommended when the genome is big, such as human and mouse and there are many motifs to examine.
}
}
\details{
This function integrates the complete pipeline for analyzing spatiotemporal
binding activities of transcription factors. It handles both simple cases, such
as one sample or one motif, and complex cases, such as multiple samples and
motifs. It imports location information of peak and motif sites and aligned
reads and outputs probabilities of each motif in each sample as false discovery
rates and inner products. When multiple samples are provided, it also outputs
clustering results of motif sites and enrichment of motifs in each cluster.
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
Clustering information is written in a file named "motif_filter_reads_cluster
.txt". Each row is a motif site and each column is a sample. The last three
columns are cluster id, motif id and motif site id.
Motif enrichment information is written in a file named "motif_cluster_
enrichment.txt". Each row represents a motif site and each column represents
either counts of motif sites, fold changes, p values and adjusted p values from
one sample.
}
\author{
Zheng Kuang
}
\examples{
#Import reads, motif sites and call peaks first
library(BayesPeak)
histonelist=vector("list",length=2)
histonelist[[1]]=data.frame(c(system.file("extdata","read1.txt",package="DynaMO",mustWork=TRUE),system.file("extdata","read2.txt",package="DynaMO",mustWork=TRUE)))
histonelist[[2]]=data.frame(c(system.file("extdata","read3.txt",package="DynaMO",mustWork=TRUE),system.file("extdata","read4.txt",package="DynaMO",mustWork=TRUE)))
motif=vector(length=5)
for(i in 1:5){
	motif[i]=system.file("extdata",paste("motif_",i,".txt",sep=""),package="DynaMO",mustWork=TRUE)
}
peak=vector("list",length=2)
for(i in 1:2){
    peak[[i]]=vector("list",length=nrow(histonelist[[i]]))
    for(j in 1:nrow(histonelist[[i]])){
        tempreads=get.reads(histonelist[[i]][j,1],150,"txt")
        tempreads1=as.data.frame(tempreads)
        tempreads2=cbind(tempreads1[,1:3],tempreads1[,5])
        colnames(tempreads2)=c("chr","start","end","strand")
        temppeak=summarize.peaks(bayespeak(tempreads2))
        peak[[i]][[j]]=GRanges(seqnames=Rle(as.character(space(temppeak))),ranges=IRanges(start=start(temppeak),end=end(temppeak)))
    }
}
#DynaMO function
DynaMO(histonelist,peak,motif,"l",2,readsmem=TRUE,150,"txt",250,20,0,0.01)
}

\keyword{ DynaMO }
