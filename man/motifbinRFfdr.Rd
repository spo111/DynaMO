\name{motifbinRFfdr}
\alias{motifbinRFfdr}
\title{
Calculate the false discovery rates of motif sites
}
\description{
Calculate the false discovery rates of motif sites as real binding sites
}
\usage{
motifbinRFfdr(histonemotifRFid, motifsegs, markernum, binnum)
}
\arguments{
  \item{histonemotifRFid}{
A list object with length equal to the number of motifs. Each element is a
vector of numbers of peaks overlapping motif sites from this motif.
}
  \item{motifsegs}{
A list object where each element is a GRanges object storing the locations of
motif sites from one motif.
}
  \item{markernum}{
The number of histone markers.
}
  \item{binnum}{
The number of bins spanning one motif site.
}
}
\details{
This function calculates the false discovery rates of motif sites as real
binding sites. It constructs a random forest model for each motif and another
random forest model for all motifs. It estimates a null distribution of votes
from background motif sites and calculates fdr-adjusted p values for each motif
sites.
}
\value{
A list object of 2 lists and the first stores the raw p values and the second
stores the fdr-adjusted p values. Each of the second-level lists is a list with
length equal to the number of motifs. Each element is a matrix of raw or fdr-
adjusted p values for each motif. Each row is a motif site and each column is
a sample.
}
\author{
Zheng Kuang
}
\examples{
#Prepare bin count
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
DynaMO(histonelist,peak,motif,"l",1,readsmem=TRUE,150,"txt",250,20,0,0.01)
#Prepare motif site segments
motifseg=vector("list",length=5)
get.motifmulti=function(x){
    result=get.motif(motif[x],150)
    return(result)
}
motifseg=mclapply(1:5,get.motifmulti,mc.cores=1)
#Prepare motif id
reads=vector("list",length=2)
for(i in 1:2){
    reads[[i]]=vector("list",length=2)
    for(j in 1:2){
        reads[[i]][[j]]=get.reads(histonelist[[i]][j,1],150,"txt")
    }
}
motifpeakoverlap=vector("list",length=2)
for(markcount in 1:2){
    motifpeakoverlap[[markcount]]=vector("list",length=length(motifseg))
    for(i in 1:length(motifseg)){
        motifpeakoverlap[[markcount]][[i]]=countOverlaps(motifseg[[i]],peak[[markcount]][[1]])
            for(j in 2:2){
                motifpeakoverlap[[markcount]][[i]]=cbind(motifpeakoverlap[[markcount]][[i]],countOverlaps(motifseg[[i]],peak[[markcount]][[j]]))
            }
    }
}
motifid0=vector("list",length=5)
motifsegs=motifseg
for(i in 1:5){
    motifid0[[i]]=1:length(motifsegs[[i]])
}
histonemotifsegreads=vector("list",length=2)
for(markcount in 1:2){
    histonemotifsegreads[[markcount]]=countmotifsegreads(reads[[markcount]],motifsegs,150,1,"txt",motifid0)
}
motifid1=motifRFid(motifsegs,histonemotifsegreads,motifpeakoverlap,motifid0)
#motifbinRFfdr
motifrandomforestfdr=motifbinRFfdr(motifid1,motifsegs,2,25)
#innerprocal
RFinnerprod=innerprocal(motifrandomforestfdr[[2]],0.01,histonemotifsegreads,1:2,25)
}

\keyword{ motifbinRFfdr }
