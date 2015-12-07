
##Abstract

Biological processes are usually associated with genome-wide remodeling of transcription driven by transcription factors (TFs). Identifying key TFs and their spatiotemporal binding patterns are indispensable to understanding how dynamic processes are programmed. Here we present a R package, dynamic motif occupancy (DynaMO), which exploits random forest modeling and clustering based enrichement analysis. DynaMO exploits TF motifs and dynamic ChIP-seq data of chromatin surrogates such as histone modifications to infer important TFs and their spatiotemporal binding in biological processes. This vignette describes how to analyze TF binding patterns using DynaMO.

##Introduction

Transcription factors (TFs) are a class of proteins that bind to functional regulatory DNA sequences, and regulate the expression of target genes. They play critical roles in ensuring the accuracy and specificity of transcription, not only under homeostatic conditions but also in various dynamic processes, such as development, differentiation and stress response. Therefore, examining the dynamics and context of TF binding and activity is important for understanding their physiological functions. Chromatin immunoprecipitation followed by deep sequencing (ChIP-seq) has been widely used to examine genome-wide binding of TFs. However, with no prior information, screening hundreds of TFs across multiple conditions is laborious and costly, and often impossible due to limitations in materials and reagents.

TFs usually recognize specific patterns of DNA sequences, known as TF binding motifs, which are characterized by consensus sequences or position-specific frequency matrices (PSFMs). Multiple computational methods have been developed to predict TF binding sites by mapping these motifs to genomes. However, methods based on motif sequences only often have high false discovery rates due to genome complexity and the relative simplicity and degeneracy of motifs. This has been attacked by incorporating genomic data on chromatin states and structures, including histone occupancy, histone modifications and DNase sensitivity. Extending this to time courses is likely to add substantial new information and reduce noise.

Here we introduce a method, called dynamic motif occupancy (DynaMO), which extends this strategy to dynamic processes. DynaMO first predicts motif binding sites of individual TFs and then clusters predicted binding sites of all TFs to characterize temporal patterns. DynaMO is distinguished from previous methods by three aspects. First, DynaMO uses a random forest model to predict where TFs bind, which is suggested to be simple and robust compared to other algorithms. Second, DynaMO predicts when TFs bind by incorporating temporal patterns of epigenomics data, such as histone modification. Third, DynaMO identifies candidate TFs that are potentially important in a dynamic process by automatic clustering and enrichment analysis. This could be very useful to screen hundreds of TFs in a poorly-studied or brand new system. DynaMO is customized to exploit time-course ChIP-seq data to better characterize the dynamics of TF network in a temporal process. 

## A simple workflow

Two types of inputs are required: (1) PSFMs or consensus sequences which represent specific sequences TFs recognize; (2) alignments of ChIP-seq reads. We used ChIP-seq of histone modification as an example but the types of epigenomic data can be extended to histone occupancy, DNase sensitivity and so on. TF motif matrices or consensus sequences are mapped to the corresponding genome using the CisGenome software package to obtain TF motif sites.

Load the package as follows:
```{r, eval=FALSE}
library(DynaMO,quietly=T)
```

The paths of read alignment files from multiple samples of one chromatin surrogate are stored in a data.frame object. Users can either read a local file storing the paths or directly input the paths in R. Paths for different epigenetic surrogates are combined in a list. The format of read alignment can be bam or txt. A txt file contains three columns separated by tab. The first is the chromosome id, the second is the coordinate and the third is the strand ("+" or "-"). In this example, we directly input the paths in R.
```{r, eval=FALSE}
histonelist<-vector("list",length=2)
histonelist[[1]]<-data.frame(c(system.file("extdata","read1.txt",package="DynaMO",mustWork=TRUE),
system.file("extdata","read2.txt",package="DynaMO",mustWork=TRUE)))
histonelist[[2]]<-data.frame(c(system.file("extdata","read3.txt",package="DynaMO",mustWork=TRUE),
system.file("extdata","read4.txt",package="DynaMO",mustWork=TRUE)))
```

The paths of motif sites files are stored in a vector object. Users can directly input the paths in R. Motif sites are stored in txt files. A txt file contains four columns separated by tab. The first is the motif note, the second is the chromosome id, the third is the left coordinate and the fourth is the right coordinate.In this example, we input 5 example motifs in R.
```{r, eval=FALSE}
motif<-vector(length=5)
for(i in 1:5){
	motif[i]<-system.file("extdata",paste("motif_",i,".txt",sep=""),package="DynaMO",mustWork=T)
}
```

Next, we need to detect peaks from each sample of aligned reads. Users can either call peaks elsewhere and input the peak regions in R or directly call peaks in R. In this example, we use BayesPeak to call peaks.
```{r, eval=FALSE}
library(BayesPeak)
peak<-vector("list",length=2)
for(i in 1:2){
    peak[[i]]<-vector("list",length=nrow(histonelist[[i]]))
    for(j in 1:nrow(histonelist[[i]])){
        tempreads<-get.reads(histonelist[[i]][j,1],150,"txt")
        tempreads1<-as.data.frame(tempreads)
        tempreads2<-cbind(tempreads1[,1:3],tempreads1[,5])
        colnames(tempreads2)<-c("chr","start","end","strand")
        temppeak<-summarize.peaks(bayespeak(tempreads2))
        peak[[i]][[j]]<-GRanges(seqnames=Rle(as.character(space(temppeak))),
                      ranges=IRanges(start=start(temppeak),end=end(temppeak)))
    }
}
```

With all input data, we can run DynaMO directly.
```{r, eval=FALSE}
DynaMO(readlist=histonelist,peak=peak,motif=motif,mode="l",core=2,readsmem=T,readlen=150,
readformat="txt",motiflen=300,motifbin=20,cluster=0,fdrcut=0.01)
```

DynaMO will automatically write the results in the current directory. Normalized read counts at motif sites are written in files named "motif_[motif index]_readcount.txt". Each row is a motif site and named by the index of motif sites. Each column is a sample. False discovery rates are written in files named "motif_[motif index]_fdr.txt". Rows are motif sites with the same order as readcount files and columns represent samples. Inner products are written in files named "motif_[motif index]_innerprod.txt". Rows are motif sites with the same order as readcount files and columns represent samples. Clustering information is written in a file named "motif_filter_reads_cluster.txt". Each row is a motif site and each column is a sample. The last three columns are cluster id, motif id and motif site id. Motif enrichment information is written in a file named "motif_cluster_enrichment.txt". Each row represents a motif site and each column represents either counts of motif sites, fold changes, p values and adjusted p values from one sample.

Some important parameters to adjust:

mode: When one or a few motifs are examined, "l" is recommended because DynaMO will evaluate every motif sites. When hundreds of motifs are examined in a small genome, such as yeast genome, "l" is also preferred because each motif only has a few thousands of motif sites from our experiences. However, when hundreds of motifs are examined in a big genome, such as human genome, "s" is recommended because each motif has hundreds of thousands of motif sites and examing motif sites close to chromatin surrogate peaks will save a lot of time.

motiflen: The length to extend from the center of motif sites depends on the type of chromatin surrogates, the method to detect surrogate peaks and the genome to study. The goal is that by extension, motif sites in the chromatin open regions will overlap chromatin surrogate peaks. Take histone modifications in yeast for example, the chromatin open region is usually 200-300 bp wide and close to a histone modification peak. When we extend motif sites by 300 bp on both sides, we will gurantee that motif sites will overlap histone modification peak.

motifbin: The length of motif bins depends on motiflen. When we extend motif sites by 200-300 bp, motifbin can be set to be 20. When motiflen is 500-1000, we can set motifbin to be 50, which will save a lot of computation.

fdrcut: fdrcut is used to select positive motif sites for inner product calculation and clustering. Default is 0.01. Users can lower the fdrcut to have more stringent motif sites or increase the fdrcut to allow more motif sites. A recommended way is to examine the distribution of FDRs and choose a proper cutoff.

## A step-by-step workflow

Alternatively, the DynaMO analysis can also be divided into separate steps and users can easily examine intermediate values and change parameters at each step. We now go into the workflow in more depth.

### Data preprocessing

Same as the previous session, we prepare the data input and get three objects, histonelist, motif and peak. Next, we input some basic parameters. "markernum" is the number of chromatin surrogates in the study. "time" is the number of sequencing samples for each surrogates. Motif sites are extended by the length of "motiflen" from the centers and each motif site is divided into consecutive bins. "motifbin" is the length of base pairs for each bin. 

```{r, eval=FALSE}
markernum<-length(histonelist)
time<-nrow(histonelist[[1]])
motifnum<-length(motif)
readlen<-150
motiflen<-300
core<-2
motifbin<-20
format<-"txt"
motifname<-1:5
```

Then we read aligned reads and store the locations in GRanges objects in R. If the memory is small, users can skip this step and just use the read list object.

```{r, eval=FALSE}
reads<-vector("list",length=markernum)
for(i in 1:markernum){
    reads[[i]]<-vector("list",length=time)
    for(j in 1:time){
        reads[[i]][[j]]<-get.reads(histonelist[[i]][j,1],readlen,"txt")
    }
}
```

Next we read motif sites and store the locations in GRanges objects in R. We also calculated the number of bins in each motif site.

```{r, eval=FALSE}
motifseg<-vector("list",length=motifnum)
get.motifmulti<-function(x){
    result=get.motif(motif[x],motiflen)
    return(result)
}
motifseg<-mclapply(1:motifnum,get.motifmulti,mc.cores=core)
motifbinnum<-(end(motifseg[[1]])[1]-start(motifseg[[1]])[1]+1)/motifbin
```

### Characterization of motif sites

TFs usually bind to open chromatin regions. Therefore, a simple way to categorize motif sites is to examine whether a motif site is overlapping or close to a peak of chromatin surrogates. So we first identify motif sites which overlap peaks.

```{r, eval=FALSE}
motifpeakoverlap<-vector("list",length=markernum)
for(markcount in 1:markernum){
    motifpeakoverlap[[markcount]]<-vector("list",length=length(motifseg))
    for(i in 1:length(motifseg)){
        motifpeakoverlap[[markcount]][[i]]<-countOverlaps(motifseg[[i]],peak[[markcount]][[1]])
        if(time>1){
            for(j in 2:time){
                motifpeakoverlap[[markcount]][[i]]<-cbind(motifpeakoverlap[[markcount]][[i]],
                countOverlaps(motifseg[[i]],peak[[markcount]][[j]]))
            } 
        }
    }
}
```

Sometimes we need to evaluate hundreds of motifs and each motif has hundreds of thousands of motif sites. In this case, evaluating motif sites which are close to chromatin open regions will save a lot of time and still keep high sensitivity.

```{r, eval=FALSE}
motifsegs<-vector("list",length=motifnum)
motifid0<-vector("list",length=motifnum)
for(i in 1:motifnum){
      temp<-rowSums(motifpeakoverlap[[1]][[i]])
      if(markernum>1){
          for(motifcount in 2:markernum){
              temp<-temp+rowSums(motifpeakoverlap[[motifcount]][[i]])
            }
          }
      motifid0[[i]]<-which(temp>0)
      motifsegs[[i]]<-motifseg[[i]][motifid0[[i]]]
}
```

Otherwise, we will examine all motif sites.

```{r, eval=FALSE}
motifsegs<-vector("list",length=motifnum)
motifid0<-vector("list",length=motifnum)
motifsegs<-motifseg
for(i in 1:motifnum){
   motifid0[[i]]<-1:length(motifsegs[[i]])
}
```

Another important feature of motif sites is the signal intensity of chromatin surrogates at motif sies. 

```{r, eval=FALSE}
histonemotifsegreads<-vector("list",length=markernum)
for(markcount in 1:markernum){
  histonemotifsegreads[[markcount]]<-countmotifsegreads(reads[[markcount]],motifsegs,readlen,
  core,format,motifid0)
}
```

### Select motif sites for random forest modeling

Based on the characterization of motif sites, we next select positive, negative and background motif sites for training random forest models. Positive and negative motif sites are used for training random forest models and background motif sites are used for estimating the null distribution.

```{r, eval=FALSE}
motifid1<-motifRFid(motifsegs,histonemotifsegreads,motifpeakoverlap,motifid0)
```

### Characterizing spatial distribution of chromatin surrogates across motif sites

We first divide each motif site into consecutive bins and store the locations in GRanges objects. Next, we count signals of surrogates at each bin.

```{r, eval=FALSE}
get.motifbintemp<-function(motifseg,bin=motifbin){
    temp<-get.motifbin(motifseg,bin)
    return(temp)
}
motifbintemp<-mclapply(motifsegs,get.motifbintemp,mc.cores=core)
for(markcount in 1:markernum){
    motifbincount(motifbintemp,reads[[markcount]],core,readlen,format,markcount)
}
```

Similarly, we perform the same procedure for motif sites selected for random forest modeling.

```{r, eval=FALSE}
motifsegrf<-vector("list",length=motifnum)
for(i in 1:motifnum){
    motifsegrf[[i]]<-vector("list",length=time)
    for(j in 1:time){
        motifsegrf[[i]][[j]]<-vector("list",length=3)
        for(k in 1:3){
            motifsegrf[[i]][[j]][[k]]<-motifseg[[i]][motifid1[[i]][[j]][[k]]]}
    }
}
get.motifbinv2<-function(motifseg,bin=motifbin){
      time<-length(motifseg)
      segnum<-length(motifseg[[1]])
      motifbinGR<-vector("list",length=time)
      for(i in 1:time){
        motifbinGR[[i]]<-vector("list",length=segnum)
        for(timenum in 1:segnum){
          chrtemp<-as.character(seqnames(motifseg[[i]][[timenum]]))
          starttemp<-start(motifseg[[i]][[timenum]])
          endtemp<-end(motifseg[[i]][[timenum]])
          motiflocnum<-length(chrtemp)
          motifbinnum<-(endtemp[1]-starttemp[1]+1)/bin
          startmatrix<-matrix(NA,nrow=motiflocnum,ncol=motifbinnum)
          endmatrix<-matrix(NA,nrow=motiflocnum,ncol=motifbinnum)
          startmatrix[,1]<-starttemp
          endmatrix[,motifbinnum]<-endtemp
          for(j in 2:motifbinnum){
            startmatrix[,j]<-startmatrix[,(j-1)]+bin
            endmatrix[,(motifbinnum-j+1)]<-endmatrix[,(motifbinnum-j+2)]-bin
          }
          motifbinGR[[i]][[timenum]]<-GRanges(seqnames=Rle(rep(chrtemp,each=motifbinnum)),
          ranges=IRanges(start=c(t(startmatrix)),end=c(t(endmatrix))))}}
      return(motifbinGR)
    }
motifbinrftemp<-mclapply(motifsegrf,get.motifbinv2,mc.cores=core)
for(markcount in 1:markernum){
      motifbincount(motifbinrftemp,reads[[markcount]],core,readlen,format,markcount)
}    
```

### Calculate FDRs and inner products

After we prepare all input data, we next calculate the FDRs and inner products of motif sites. An important parameter to adjust is fdrcut. Motif sites with FDR less than fdrcut are used as reference to calculate inner products. Here we use 0.01.

```{r, eval=FALSE}
motifrandomforestfdr<-motifbinRFfdr(motifid1,motifsegs,markernum,motifbinnum)
RFinnerprod<-innerprocal(motifrandomforestfdr[[2]],0.01,histonemotifsegreads,1:markernum,
motifbinnum) 
```

We can use FDR to select positive motif sites and use inner products to further rank motif sites when FDRs are same.

### Clustering of motif sites

When multiple samples are provided, we can cluster positive motif sites based on the patterns of signals of chromatin surrogates at motif sites. In the previous session, we have already calculated signal intensity of chromatin surrogates and store it in the histonemotifsegreads object. We now need to combine the signals at positive motif sites from different motifs. First, we find the motif sites with FDR less than fdrcut and combine read counts at these motif sites from all motifs in a large matrix.

```{r, eval=FALSE}
fdrcut<-0.00001
motiffdrread<-matrix(ncol=time*markernum)
for(i in 1:motifnum){
  temp<-t(scale(t(histonemotifsegreads[[1]][[i]])))
  if(markernum>1){
    for(markcount in 2:markernum){
      temp<-cbind(temp,t(scale(t(histonemotifsegreads[[markcount]][[i]]))))
    }
  }
  rownames(temp)<-paste(i,rownames(temp),sep="_")
  motifrandomforestfdr[[2]][[i]][is.na(motifrandomforestfdr[[2]][[i]])]=1
  numtest<-apply(motifrandomforestfdr[[2]][[i]],1,function(x){if(min(x)<=fdrcut){return(1)}
    else{return(0)}})
  motiffdrread<-rbind(motiffdrread,temp[numtest==1,])
}
motiffdrread<-na.exclude(motiffdrread)
```

Next, we perform Kmeans clustering. Users can either do automatic Kmeans or choose a specific value for the number of clusters. Here we use fpc to automatically generate the number of clusters.

```{r, eval=FALSE}
library(fpc)
cluster<-0
#cluster=0 means automatic clustering
if(cluster==0){
        tempread<-motiffdrread[sample(1:nrow(motiffdrread),min(c(nrow(motiffdrread),10000))),]
        tempcls<-pamk(tempread,2:20,usepam=F)
        clsnum<-tempcls$nc
        }else{
            clsnum<-cluster
        }
        motifk<-kmeans(motiffdrread,centers=clsnum,10,10)
        motiffdrread1<-cbind(motiffdrread,motifk$cluster)
        tempname<-matrix(as.numeric(unlist(strsplit(rownames(motiffdrread1),split="_"))),
                         ncol=2,byrow=T)
        motiffdrread1<-cbind(motiffdrread1,tempname)
        colnames(motiffdrread1)<-c(paste(rep(paste("marker",1:markernum,sep=""),each=time),
                                         1:time,sep="_"),
        "cluster","motifname","motiflocid")
```

The object motiffdrread1 contains read counts of chromatin surrogates at motif sites and clustering information of each motif site. The pattern of read counts at a motif site represents the binding pattern of the corresponding TF at this region.

Next, we perform fisher.exact test to find motifs with motif sites enriched in the clusters.

```{r, eval=FALSE}
motiffdrclstest=motifclstest(motiffdrread1,1:motifnum)
colnames(motiffdrclstest)=paste(rep(paste("cluster",1:clsnum)),rep(c("count","ratio","p-value","padj"),each=4))
```

Motifs enriched in any of the clusters are predicted to be important and the averaged pattern of read counts in the cluster which a motif is enriched reflects when the TF is mainly functional.

