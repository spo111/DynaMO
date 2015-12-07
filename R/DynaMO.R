
get.reads<-function(file,seqL=150,readformat="bam"){
    if(readformat=="bam"){
      reads<-readGAlignments(as.character(file))
      reads1<-reads[strand(reads)=="+"]
      reads2<-reads[strand(reads)=="-"]
      res<-GRanges(seqnames=c(seqnames(reads1),seqnames(reads2)),ranges=c(IRanges(start=start(reads1),width=seqL),IRanges(end=end(reads2),width=seqL)),strand=c(strand(reads1),strand(reads2)))
    }else{
        reads<-read.table(as.character(file),sep="\t",col.names=c("chrom","position","strand"))
        reads2<-split(reads,reads$strand)
        res<-GRanges(seqnames=Rle(c(as.character(reads2[["+"]]$chrom),as.character(reads2[["-"]]$chrom))),ranges=c(IRanges(start=reads2[["+"]]$position,width=seqL),IRanges(end=reads2[["-"]]$position,width=seqL)),strand=Rle(c("+","-"),sapply(reads2,nrow)))
    }
    return(res)
}

get.motif<-function(file,extendL=300){
    motif<-try(read.table(as.character(file),sep="\t"),silent=T)
    if(is.data.frame(motif)){
        motifcenter=as.integer(rowMeans(motif[,3:4]))
        motifGR<-GRanges(seqnames=Rle(as.character(motif[,2])),ranges=IRanges(start=(motifcenter-extendL),end=(motifcenter+extendL-1)))}
}

get.motifbin<-function(motifseg,bin=50){
    chrtemp=as.character(seqnames(motifseg))
    starttemp=start(motifseg)
    endtemp=end(motifseg)
    motiflocnum=length(chrtemp)
    motifbinnum=(endtemp[1]-starttemp[1]+1)/bin
    startmatrix=matrix(NA,nrow=motiflocnum,ncol=motifbinnum)
    endmatrix=matrix(NA,nrow=motiflocnum,ncol=motifbinnum)
    startmatrix[,1]=starttemp
    endmatrix[,motifbinnum]=endtemp
    for(i in 2:motifbinnum){
        startmatrix[,i]=startmatrix[,(i-1)]+bin
        endmatrix[,(motifbinnum-i+1)]=endmatrix[,(motifbinnum-i+2)]-bin
    }
    motifbinGR=GRanges(seqnames=Rle(rep(chrtemp,each=motifbinnum)),ranges=IRanges(start=c(t(startmatrix)),end=c(t(endmatrix))))
    return(motifbinGR)
}
countmotifsegreads<-function(read,segmentfile,seqL=150,core=1,format="bam",histonemotifRFid){
  countmotifsegreadsv1<-function(listfile,segmentfile,seqL,core,format="bam",histonemotifRFid){
    segnum=length(segmentfile)
    result<-vector("list",segnum)
    for(i in 1:segnum){
      result[[i]]<-matrix(NA,nrow=length(segmentfile[[i]]),ncol=nrow(listfile))
    }
    totalreads<-vector(length=nrow(listfile))
    for(i in 1:nrow(listfile)){
      temp<-get.reads(listfile[i,],seqL,format)
      totalreads[i]=length(temp)
      temp1=vector("list",length=core)
      if(core==1){
        temp1[[1]]=temp
        rm(temp)
      }else{
        readcount=totalreads[i]
        readseg=round(readcount/core)
        for(j in 1:(core-1)){
          temp1[[j]]=temp[((j-1)*readseg+1):(j*readseg)]
        }
        temp1[[core]]=temp[((core-1)*readseg+1):readcount]
        rm(temp)
      }
      for(j in 1:segnum){
        countreads=function(x){
          tempresult=try(countOverlaps(segmentfile[[j]],temp1[[x]]))
          return(tempresult)
        }
        tempbincount=mclapply(1:core,countreads,mc.cores=core)
        tempbincount1=tempbincount[[1]]
        if(core>1){
          for(k in 2:core){
            tempbincount1=as.numeric(tempbincount1)+as.numeric(tempbincount[[k]])
          }
        }
        rm(tempbincount)
        result[[j]][,i]=tempbincount1
      }
    }
    ratio=totalreads/mean(totalreads)
    for(j in 1:segnum){
      for(i in 1:nrow(listfile)){
        result[[j]][,i]=try(round(log2(result[[j]][,i]/ratio[i]+1),5),silent=T)
      }
      if(nrow(result[[j]]>0))
      {rownames(result[[j]])=try(histonemotifRFid[[j]],silent=T)
      }    }
    return(result)
  }
  countmotifsegreadsv2<-function(reads,segmentfile,core,histonemotifRFid){
    segnum=length(segmentfile)
    result<-vector("list",segnum)
    for(i in 1:segnum){
      result[[i]]<-matrix(NA,nrow=length(segmentfile[[i]]),ncol=length(reads))
    }
    totalreads<-sapply(reads,length)
    for(i in 1:length(reads)){
      temp1=vector("list",length=core)
      if(core==1){
        temp1[[1]]=reads[[i]]
      }else{
        readseg=round(totalreads[i]/core)
        for(j in 1:(core-1)){
          temp1[[j]]=reads[[i]][((j-1)*readseg+1):(j*readseg)]
        }
        temp1[[core]]=reads[[i]][((core-1)*readseg+1):totalreads[i]]
      }
      for(j in 1:segnum){
        countreads=function(x){
          tempresult=try(countOverlaps(segmentfile[[j]],temp1[[x]]))
          return(tempresult)
        }
        tempbincount=mclapply(1:core,countreads,mc.cores=core)
        tempbincount1=as.numeric(tempbincount[[1]])
        if(core>1){
          for(k in 2:core){
            tempbincount1=as.numeric(tempbincount1)+as.numeric(tempbincount[[k]])
          }
        }
        rm(tempbincount)
        result[[j]][,i]=tempbincount1
      }
      rm(temp1)
      gc()
    }
    ratio=totalreads/mean(totalreads)
    for(j in 1:segnum){
      for(i in 1:length(totalreads)){
        result[[j]][,i]=try(round(log2(result[[j]][,i]/ratio[i]+1),5))
      }
      if(nrow(result[[j]]>0))
      {rownames(result[[j]])=histonemotifRFid[[j]]
      }  }
    return(result)
  }
    if(is.list(read)){
        result=countmotifsegreadsv2(read,segmentfile,core,histonemotifRFid)
    }else{
        result=countmotifsegreadsv1(read,segmentfile,seqL,core,format,histonemotifRFid)
    }
    return(result)
}


motifbincount=function(motifbin,reads,core=1,seqL=150,format="bam",name){
    segnum=length(motifbin)
    if(is.list(reads)){
        time=length(reads)
        readtype=T
    }else{
        time=nrow(reads)
        readtype=F
    }
    if(is.list(motifbin[[1]])){
        motiftype=T
    }else{
        motiftype=F
    }
    for(i in 1:time){
        if(readtype){
            temp=reads[[i]]
        }else{
            temp<-get.reads(listfile[i,],seqL,format)
        }
        temp1=vector("list",length=core)
        readcount=length(temp)
        readseg=round(readcount/core)
        if(core>1){
          for(j in 1:(core-1)){
            temp1[[j]]=temp[((j-1)*readseg+1):(j*readseg)]
          }
          temp1[[core]]=temp[((core-1)*readseg+1):readcount]
        }else{
          temp1[[1]]=temp
        }
        rm(temp)
        for(j in 1:segnum){
            if(motiftype){
                for(k in 1:3){
                    if(length(motifbin[[j]][[i]][[k]])>0){
                      countreads=function(x){
                        tempresult=try(countOverlaps(motifbin[[j]][[i]][[k]],temp1[[x]]))
                        return(tempresult)
                      }
                      tempbincount=mclapply(1:core,countreads,mc.cores=core)
                      tempbincount1=tempbincount[[1]]
                      if(core>1){
                        for(l in 2:core){
                          tempbincount1=tempbincount1+tempbincount[[l]]}
                      }
                      rm(tempbincount)
                      tempbincount1=round(log2(tempbincount1+1),5)
                    }else{tempbincount1=0}
                    write.table(tempbincount1,file=paste("motifbin_",name,"_",i,"_",j,"_",k,".txt",sep=""),quote=F,sep="\t",row.names=F,col.names=F)
                }
            }else{
                if(length(motifbin[[j]])>0){
                  countreads=function(x){
                    tempresult=try(countOverlaps(motifbin[[j]],temp1[[x]]))
                    return(tempresult)
                  }
                  tempbincount=mclapply(1:core,countreads,mc.cores=core)
                  tempbincount1=tempbincount[[1]]
                  if(core>1){
                    for(k in 2:core){
                      tempbincount1=tempbincount1+tempbincount[[k]]
                    }
                  }
                  rm(tempbincount)
                  tempbincount1=round(log2(tempbincount1+1),5)
                }else{tempbincount1=0}
                write.table(tempbincount1,file=paste("motifbin_",name,"_",i,"_",j,".txt",sep=""),quote=F,sep="\t",row.names=F,col.names=F)
            }
        }
        rm(temp1)
        gc()
    }}

motifRFid=function(motifseg,histonemotifsegreads,motifpeakoverlap,motifseg500overlapid){
    segnum=length(motifseg)
    marknum=length(histonemotifsegreads)
    result=vector("list",length=segnum)
    time=ncol(histonemotifsegreads[[1]][[1]])
    for(i in 1:segnum){
        result[[i]]=vector("list",length=time)
        for(j in 1:time){
            result[[i]][[j]]=vector("list",length=3)
            tempmarkreads=vector(length=nrow(histonemotifsegreads[[1]][[i]]))
            temppeakoverlap=vector(length=length(motifseg[[i]]))
            for(markcount in 1:marknum){
                tempmarkreads=tempmarkreads+histonemotifsegreads[[markcount]][[i]][,j]
                temppeakoverlap=temppeakoverlap+motifpeakoverlap[[markcount]][[i]][,j]
            }
            temp=cbind(1:length(motifseg[[i]]),temppeakoverlap)
            if(nrow(temp)>200){
                for(markcount in marknum:1){
                    if(nrow(temp[temp[,2]>=markcount,])>=250){
                        temppeak=temp[temp[,2]>=markcount,]
                        break
                    }
                }
                if(markcount==1){
                    temppeak=temp[temp[,2]>0,]
                }
                tempnopeak=temp[temp[,2]==0,]
                peakvalue=cbind(motifseg500overlapid[[i]],tempmarkreads)
                peakvalue1=peakvalue[peakvalue[,1]%in%temppeak[,1],]
                if(is.matrix(temppeak)&is.matrix(tempnopeak)){
                    if(nrow(peakvalue1)>500){
                        peakvalue1=peakvalue1[order(peakvalue1[,2],decreasing=T),]
                        if(nrow(tempnopeak)>=250){
                            temppeakid=sample(peakvalue1[1:round(nrow(peakvalue1)/2),1],250)
                            tempnopeakid=sample(tempnopeak[,1],250)
                        }else{
                            temppeakid=sample(peakvalue1[1:round(nrow(peakvalue1)/2),1],(500-nrow(tempnopeak)))
                            tempnopeakid=tempnopeak[,1]
                        }
                    }else if(nrow(peakvalue1)>250){
                        if(nrow(tempnopeak)>=250){
                            temppeakid=sample(peakvalue1[,1],250)
                            tempnopeakid=sample(tempnopeak[,1],250)
                        }else{
                            temppeakid=sample(peakvalue1[,1],min(500-nrow(tempnopeak),nrow(peakvalue1)))
                            tempnopeakid=tempnopeak[,1]
                        }
                    }else{
                        temppeakid=peakvalue1[,1]
                        tempnopeakid=sample(tempnopeak[,1],min(500-nrow(peakvalue1),nrow(tempnopeak)))
                    }
                    result[[i]][[j]][[1]]=temppeakid
                    result[[i]][[j]][[2]]=tempnopeakid
                    result[[i]][[j]][[3]]=sample(tempnopeak[,1],500,replace=T)
                }
            }
        }
    }
    return(result)
}

motifbinRFfdr=function(histonemotifRFid,motifsegs,markernum,binnum){
    segnum=length(histonemotifRFid)
    time=length(histonemotifRFid[[1]])
    result=vector("list",length=segnum)
    for(i in 1:segnum){
        result[[i]]=matrix(nrow=length(motifsegs[[i]]),ncol=time)
    }
    result1=result
    for(i in 1:time){
        tempbincounttotal=vector("list",length=3)
        temprfval=vector("list",length=segnum)
        for(j in 1:segnum){
            temprfval[[j]]=vector("list",length=4)
            tempbincount1=vector("list",length=markernum)
            for(markcount in 1:markernum){
                tempbincount=read.table(paste("motifbin_",markcount,"_",i,"_",j,".txt",sep=""),sep="\t")
                tempbincount1[[markcount]]=tempbincount[,1]
            }
            if(try(length(histonemotifRFid[[j]][[i]][[1]]))>0){
                temppeakid=read.table(paste("motifbin_",1,"_",i,"_",j,"_",1,".txt",sep=""),sep="\t")
                tempnopeakid=read.table(paste("motifbin_",1,"_",i,"_",j,"_",2,".txt",sep=""),sep="\t")
                tempbgid=read.table(paste("motifbin_",1,"_",i,"_",j,"_",3,".txt",sep=""),sep="\t")
                tempoutpeak=matrix(temppeakid[,1],byrow=T,ncol=binnum)
                tempoutnopeak=matrix(tempnopeakid[,1],byrow=T,ncol=binnum)
                bgtempbincount=matrix(tempbgid[,1],byrow=T,ncol=binnum)
                tempbincountset=matrix(tempbincount1[[1]],ncol=binnum,byrow=T)
                if(markernum>1){
                    for(markcount in 2:markernum){
                        temppeakid=read.table(paste("motifbin_",markcount,"_",i,"_",j,"_",1,".txt",sep=""),sep="\t")
                        tempnopeakid=read.table(paste("motifbin_",markcount,"_",i,"_",j,"_",2,".txt",sep=""),sep="\t")
                        tempbgid=read.table(paste("motifbin_",markcount,"_",i,"_",j,"_",3,".txt",sep=""),sep="\t")
                        tempoutpeak0=matrix(temppeakid[,1],byrow=T,ncol=binnum)
                        tempoutnopeak0=matrix(tempnopeakid[,1],byrow=T,ncol=binnum)
                        bgtempbincount0=matrix(tempbgid[,1],byrow=T,ncol=binnum)
                        tempoutpeak=cbind(tempoutpeak,tempoutpeak0)
                        tempoutnopeak=cbind(tempoutnopeak,tempoutnopeak0)
                        bgtempbincount=cbind(bgtempbincount, bgtempbincount0)
                        tempbincountset=cbind(tempbincountset,matrix(tempbincount1[[markcount]],ncol=binnum,byrow=T))
                    }
                }
                tempbincounttotal[[1]]=rbind(tempbincounttotal[[1]],tempoutpeak)
                tempbincounttotal[[2]]=rbind(tempbincounttotal[[2]],tempoutnopeak)
                tempbincounttotal[[3]]=rbind(tempbincounttotal[[3]],bgtempbincount)
                trainset=rbind(tempoutpeak,tempoutnopeak)
                trainresponse=c(rep(1,nrow(tempoutpeak)),rep(0,nrow(tempoutnopeak)))
                motifbinrfmodel=randomForest(x=trainset,y=factor(trainresponse))
                predictmotifbinrf=predict(motifbinrfmodel,newdata=tempbincountset,type="vote")
                temprfval[[j]][[1]]=predictmotifbinrf[,2]
                predictmotifbinbg=predict(motifbinrfmodel,newdata=bgtempbincount,type="vote")
                temprfval[[j]][[2]]=predictmotifbinbg[,2]
            }
        }
        temprow=min(10000,nrow(tempbincounttotal[[1]]),nrow(tempbincounttotal[[2]]))
        trainset=rbind(tempbincounttotal[[1]][sample(nrow(tempbincounttotal[[1]]),temprow),],tempbincounttotal[[2]][sample(nrow(tempbincounttotal[[2]]),temprow),])
        trainresponse=c(rep(1,temprow),rep(0,temprow))
        motifbinrfmodeltotal=randomForest(x=trainset,y=factor(trainresponse))
        for(j in 1:segnum){
            tempbincount1=vector("list",length=markernum)
            for(markcount in 1:markernum){
                tempbincount=read.table(paste("motifbin_",markcount,"_",i,"_",j,".txt",sep=""),sep="\t")
                tempbincount1[[markcount]]=tempbincount[,1]
            }
            if(try(length(histonemotifRFid[[j]][[i]][[1]]))>0){
                tempbgid=read.table(paste("motifbin_",1,"_",i,"_",j,"_",3,".txt",sep=""),sep="\t")
                bgtempbincount=matrix(tempbgid[,1],byrow=T,ncol=binnum)
                tempbincountset=matrix(tempbincount1[[1]],ncol=binnum,byrow=T)
                if(markernum>1){
                    for(markcount in 2:markernum){
                        tempbgid=read.table(paste("motifbin_",markcount,"_",i,"_",j,"_",3,".txt",sep=""),sep="\t")
                        bgtempbincount0=matrix(tempbgid[,1],byrow=T,ncol=binnum)
                        bgtempbincount=cbind(bgtempbincount, bgtempbincount0)
                        tempbincountset=cbind(tempbincountset,matrix(tempbincount1[[markcount]],ncol=binnum,byrow=T))
                    }
                }
                predictmotifbinrf=predict(motifbinrfmodeltotal,newdata=tempbincountset,type="vote")
                temprfval[[j]][[3]]=predictmotifbinrf[,2]
                predictmotifbinbg=predict(motifbinrfmodeltotal,newdata=bgtempbincount,type="vote")
                temprfval[[j]][[4]]=predictmotifbinbg[,2]
                motifbinbgecdf=ecdf((temprfval[[j]][[2]]+temprfval[[j]][[4]])/2)
                predictmotifbinrfpval=1-motifbinbgecdf((temprfval[[j]][[1]]+temprfval[[j]][[3]])/2)
                predictmotifbinrfpadj=p.adjust(predictmotifbinrfpval,method="fdr")
                result[[j]][,i]=predictmotifbinrfpval
                result1[[j]][,i]=predictmotifbinrfpadj
            }
        }
    }
    return(list(result,result1))
}

innerprocal=function(RFfdr,fdrcut=0.01,RFreads,name,binnum){
    segnum=length(RFfdr)
    time=ncol(RFfdr[[1]])
    marknum=length(RFreads)
    result=vector("list",length=segnum)
    for(i in 1:segnum){
        result[[i]]=matrix(nrow=nrow(RFfdr[[i]]),ncol=time)
    }
    for(i in 1:time){
        tempbincounttotal=vector(length=binnum*length(name))
        for(j in 1:segnum){
            tempbincount1=matrix(nrow=nrow(RFfdr[[j]]))
            for(namecount in 1:length(name)){
                tempbincount=read.table(paste("motifbin_",namecount,"_",i,"_",j,".txt",sep=""),sep="\t")
                tempbincount1=cbind(tempbincount1,matrix(tempbincount[,1],ncol=binnum,byrow=T))
            }
            tempbincountset=tempbincount1[,2:ncol(tempbincount1)]
            tempfdr=RFfdr[[j]][,i]
            tempmarkreads=vector(length=nrow(RFfdr[[j]]))
            for(markcount in 1:marknum){
                tempmarkreads=tempmarkreads+RFreads[[markcount]][[j]][,i]
            }
            tempbincountset=cbind(tempbincountset,tempmarkreads)
            tempbincountsetf=tempbincountset[tempfdr<fdrcut,]
            tempbincountsetf= tempbincountsetf[order(tempbincountsetf[,ncol(tempbincountsetf)],decreasing=T),]
            tempbincountsetf=tempbincountsetf[,1:(ncol(tempbincountsetf)-1)]
            if(nrow(tempbincountsetf)>100){
                tempbincounttotal=rbind(tempbincounttotal, tempbincountsetf[sample(1:100,50),])
            }else if(nrow(tempbincountsetf)>50){
                tempbincounttotal=rbind(tempbincounttotal, tempbincountsetf[sample(1:nrow(tempbincountsetf),50),])
            }else{
                tempbincounttotal=rbind(tempbincounttotal, tempbincountsetf[,])
            }
        }
        tempbincounttotal=tempbincounttotal[2:nrow(tempbincounttotal),]
        tempbincounttotal=na.exclude(tempbincounttotal)
        tempbincounttotal[is.infinite(tempbincounttotal)]=0
        tempbinmean=colMeans(tempbincounttotal)
        for(j in 1:segnum){
            tempbincount1=matrix(nrow=nrow(RFfdr[[j]]))
            for(namecount in 1:length(name)){
                tempbincount=read.table(paste("motifbin_",namecount,"_",i,"_",j,".txt",sep=""),sep="\t")
                tempbincount1=cbind(tempbincount1,matrix(tempbincount[,1],ncol=binnum,byrow=T))
            }
            tempbincountset=tempbincount1[,2:ncol(tempbincount1)]
            for(k in 1:nrow(tempbincountset)){
                result[[j]][k,i]=sum(tempbinmean*tempbincountset[k,])

            }
        }
    }
    return(result)
}

motifclstest<-function(motifcls,motifname){
    motiflen=length(motifname)
    clsnum=length(table(motifcls[,ncol(motifcls)-2]))
    motifclscount=matrix(nrow=motiflen,ncol=3*clsnum)
    motifclstotal=table(motifcls[,ncol(motifcls)-2])
    motifclssum=nrow(motifcls)
    for(i in 1:motiflen){
        if(length(motifcls[motifcls[,ncol(motifcls)-1]==motifname[i],1])>0){
            motifclstemp=motifcls[motifcls[,ncol(motifcls)-1]==motifname[i],]
            motifclstemptotal=nrow(motifclstemp)
            for(j in 1:clsnum){
                motifclscount[i,j]=length(motifclstemp[motifclstemp[,ncol(motifcls)-2]==j,1])
                motifclscount[i,(j+2*clsnum)]=fisher.test(matrix(c(motifclscount[i,j],motifclstemptotal-motifclscount[i,j],motifclstotal[j],motifclssum-motifclstotal[j]),nrow=2),alternative="greater")$p.value
                motifclscount[i,(j+clsnum)]=(motifclscount[i,j]/motifclstemptotal)/(motifclstotal[j]/motifclssum)
            }
        }
    }
    rownames(motifclscount)=as.character(motifname)
    temp=matrix(p.adjust(as.vector(motifclscount[,(2*clsnum+1):(3*clsnum)]),"bonferroni"),ncol=clsnum)
    motifclscount=cbind(motifclscount,temp)
    return(motifclscount)
}



DynaMO<-function(readlist,peak,motif,mode="l",core=1,readsmem=TRUE,readlen=150,readformat="bam",motiflen=300,motifbin=20,cluster=0,fdrcut=0.01,batch=T){
  DynaMOs<-function(reads,peak,motifseg,core,readsmem,markernum,time,mode,format,motifbin,readlen,motifname){
    motifnum=length(motifseg)
    motifbinnum=(end(motifseg[[1]])[1]-start(motifseg[[1]])[1]+1)/motifbin
    motifpeakoverlap=vector("list",length=markernum)
    for(markcount in 1:markernum){
      motifpeakoverlap[[markcount]]=vector("list",length=length(motifseg))
      for(i in 1:length(motifseg)){
        motifpeakoverlap[[markcount]][[i]]=data.matrix(countOverlaps(motifseg[[i]],peak[[markcount]][[1]]))
        if(time>1){
          for(j in 2:time){
            motifpeakoverlap[[markcount]][[i]]=cbind(motifpeakoverlap[[markcount]][[i]],countOverlaps(motifseg[[i]],peak[[markcount]][[j]]))
          }
        }
      }
    }
    print("motifpeakoverlap")
    motifsegs=vector("list",length=motifnum)
    motifid0=vector("list",length=motifnum)
    if(mode=="l"){
      motifsegs=motifseg
      for(i in 1:motifnum){
        motifid0[[i]]=1:length(motifsegs[[i]])
      }
    }else{
      for(i in 1:motifnum){
        temp=rowSums(motifpeakoverlap[[1]][[i]])
        if(markernum>1){
          for(motifcount in 2:markernum){
            temp=temp+rowSums(motifpeakoverlap[[motifcount]][[i]])
          }
        }
        motifid0[[i]]=which(temp>0)
        motifsegs[[i]]=motifseg[[i]][motifid0[[i]]]
      }
    }
    histonemotifsegreads=vector("list",length=markernum)
    for(markcount in 1:markernum){
      histonemotifsegreads[[markcount]]=countmotifsegreads(reads[[markcount]],motifsegs,readlen,core,format,motifid0)
    }
    print("motifsegreadcount")
    motifid1=motifRFid(motifseg,histonemotifsegreads,motifpeakoverlap,motifid0)
    print("motifRFid")
    get.motifbintemp<-function(motifseg,bin=motifbin){
      temp=get.motifbin(motifseg,bin)
      return(temp)
    }
    motifbintemp=mclapply(motifsegs,get.motifbintemp,mc.cores=core)
    print("motifsegbin")
    for(markcount in 1:markernum){
      motifbincount(motifbintemp,reads[[markcount]],core,readlen,format,markcount)
    }
    print("motifbinreadcount")
    rm(motifbintemp)
    motifsegrf=vector("list",length=motifnum)
    for(i in 1:motifnum){
      motifsegrf[[i]]=vector("list",length=time)
      for(j in 1:time){
        motifsegrf[[i]][[j]]=vector("list",length=3)
        for(k in 1:3){
          motifsegrf[[i]][[j]][[k]]=motifseg[[i]][motifid1[[i]][[j]][[k]]]}
      }
    }
    get.motifbinv2<-function(motifseg,bin=motifbin){
      time=length(motifseg)
      segnum=length(motifseg[[1]])
      motifbinGR=vector("list",length=time)
      for(i in 1:time){
        motifbinGR[[i]]=vector("list",length=segnum)
        for(timenum in 1:segnum){
          if(length(motifseg[[i]][[timenum]])>0){
            chrtemp=as.character(seqnames(motifseg[[i]][[timenum]]))
            starttemp=start(motifseg[[i]][[timenum]])
            endtemp=end(motifseg[[i]][[timenum]])
            motiflocnum=length(chrtemp)
            motifbinnum=(endtemp[1]-starttemp[1]+1)/bin
            startmatrix=matrix(NA,nrow=motiflocnum,ncol=motifbinnum)
            endmatrix=matrix(NA,nrow=motiflocnum,ncol=motifbinnum)
            startmatrix[,1]=starttemp
            endmatrix[,motifbinnum]=endtemp
            for(j in 2:motifbinnum){
              startmatrix[,j]=startmatrix[,(j-1)]+bin
              endmatrix[,(motifbinnum-j+1)]=endmatrix[,(motifbinnum-j+2)]-bin
            }
            motifbinGR[[i]][[timenum]]=GRanges(seqnames=Rle(rep(chrtemp,each=motifbinnum)),ranges=IRanges(start=c(t(startmatrix)),end=c(t(endmatrix))))}
        }
      }
      return(motifbinGR)
    }
    motifbinrftemp=mclapply(motifsegrf,get.motifbinv2,mc.cores=core)
    print("motifsegbinrf")
    for(markcount in 1:markernum){
      motifbincount(motifbinrftemp,reads[[markcount]],core,readlen,format,markcount)
    }
    print("motifbinreadcountrf")
    motifrandomforestfdr=motifbinRFfdr(motifid1,motifsegs,markernum,motifbinnum)
    print("motifrandomforestfdr")
    for(i in 1:motifnum){
      motifrandomforestfdr[[2]][[i]][is.na(motifrandomforestfdr[[2]][[i]])]=1
      write.table(motifrandomforestfdr[[2]][[i]],paste("motif",motifname[i],"fdr.txt",sep="_"),quote=F,sep="\t",row.names=F,col.names=F)
    }
    RFinnerprod=innerprocal(motifrandomforestfdr[[2]],0.01,histonemotifsegreads,1:markernum,motifbinnum)
    print("RFinnerprod")
    for(i in 1:motifnum){
      write.table(RFinnerprod[[i]],paste("motif",motifname[i],"innerprod.txt",sep="_"),quote=F,sep="\t",row.names=F,col.names=F)
    }
    for(i in 1:motifnum){
      temp=t(scale(t(histonemotifsegreads[[1]][[i]])))
      if(markernum>1){
        for(markcount in 2:markernum){
          temp=cbind(temp,t(scale(t(histonemotifsegreads[[markcount]][[i]]))))
        }
      }
      write.table(round(temp,5),paste("motif_",motifname[i],"_readcount.txt",sep=""),sep="\t",row.names=T,col.names=F,quote=F)
    }
    print("motifreadwrite")
  }
    markernum=length(readlist)
    time=nrow(readlist[[1]])
    motifnum=length(motif)
    motifname=1:motifnum
    if(readsmem){
        reads=vector("list",length=markernum)
        for(i in 1:markernum){
            reads[[i]]=vector("list",length=time)
            for(j in 1:time){
                reads[[i]][[j]]=get.reads(readlist[[i]][j,1],readlen,readformat)
            }
        }
    }else{
        reads=readlist
    }
    motifseg=vector("list",length=motifnum)
    get.motifmulti=function(x){
        result=get.motif(motif[x],motiflen)
        return(result)
    }
    motifseg=mclapply(1:motifnum,get.motifmulti,mc.cores=core)
    if(motifnum>15&batch){
        motifbatch=10
        temp=motifnum%%motifbatch
        if(temp>5){
            motifround=motifnum%/%motifbatch+1
        }else{
            motifround=motifnum%/%motifbatch
        }
        for(i in 1:(motifround-1)){
            motifsegtemp=motifseg[((i-1)*motifbatch+1):(i*motifbatch)]
            motifnametemp=((i-1)*motifbatch+1):(i*motifbatch)
            DynaMOs(reads,peak,motifsegtemp,core,readsmem,markernum,time,mode,readformat,motifbin,readlen,motifnametemp)
        }
        motifsegtemp=motifseg[((motifround-1)*motifbatch+1):motifnum]
        motifnametemp=((motifround-1)*motifbatch+1):motifnum
        DynaMOs(reads,peak,motifsegtemp,core,readsmem,markernum,time,mode,readformat,motifbin,readlen,motifnametemp)
    }else{
        DynaMOs(reads,peak,motifseg,core,readsmem,markernum,time,mode,readformat,motifbin,readlen,motifname)
    }
    if(time>1){
        motiffdr=read.table("motif_1_fdr.txt",sep="\t")
        numtest=apply(motiffdr,1,function(x){if(min(x)<=fdrcut){return(1)}else{return(0)}})
        temp=try(read.table("motif_1_readcount.txt",sep="\t"))
        temp1=try(temp[,2:ncol(temp)])
        rownames(temp1)=paste(1,temp[,1],sep="_")
        motiffdrread=temp1[numtest==1,]
        if(motifnum>1){
            for(i in 2:motifnum){
                motiffdr=read.table(paste("motif",motifname[i],"fdr.txt",sep="_"),sep="\t")
                numtest=apply(motiffdr,1,function(x){if(min(x)<=fdrcut){return(1)}else{return(0)}})
                temp=try(read.table(paste("motif",i,"readcount.txt",sep="_"),sep="\t"))
                temp1=try(temp[,2:ncol(temp)])
                rownames(temp1)=paste(i,temp[,1],sep="_")
                motiffdrread=rbind(motiffdrread,temp1[numtest==1,])
            }
        }
        motiffdrread=na.exclude(motiffdrread)
        if(cluster==0){
            tempread=motiffdrread[sample(1:nrow(motiffdrread),min(c(nrow(motiffdrread),10000))),]
            tempcls=pamk(tempread,2:20,usepam=F)
            clsnum=tempcls$nc
        }else{
            clsnum=cluster
        }
        motifk=kmeans(motiffdrread,centers=clsnum,20,20)
        motiffdrread1=cbind(motiffdrread,motifk$cluster)
        tempname=matrix(as.numeric(unlist(strsplit(rownames(motiffdrread1),split="_"))),ncol=2,byrow=TRUE)
        motiffdrread1=cbind(motiffdrread1,tempname)
        colnames(motiffdrread1)=c(paste(rep(paste("marker",1:markernum,sep=""),each=time),1:time,sep="_"),"cluster","motifname","motiflocid")
        write.table(motiffdrread1,"motif_filter_reads_cluster.txt",sep="\t",row.names=F,quote=F)
        print("clustering finished")
        motiffdrclstest=motifclstest(motiffdrread1,1:motifnum)
        colnames(motiffdrclstest)=paste(rep(paste("cluster",1:clsnum)),rep(c("count","ratio","p-value","padj"),each=clsnum))
        write.table(motiffdrclstest,"motif_cluster_enrichment.txt",sep="\t",quote=FALSE)
    }
}
