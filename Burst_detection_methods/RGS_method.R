#Implement RGS method
RGS.method<-function(spike.train, thresh=0.01) {
  isi.list<-list()
  isi.list[[1]]<-data.frame(isi = c(spike.train[1], diff(spike.train)))
  if (dim(isi.list[[1]])[1]<3){
    return(NA)
  }
  bursts<-f.BPsummary(data=isi.list, Pthresh=thresh)$burst[[1]][,c(4,7:8)]
  if (dim(bursts)[1]<1){
      return(NA)
  }
  zero.bursts<-which(bursts[,2]==0) 
  if(any(zero.bursts)){
    bursts<-bursts[-zero.bursts,]
    if (dim(bursts)[1]<1){
      return(NA)
    }
  }
  combine.clusters<-function(clust.num){
  clust<-which(bursts$clusid==clust.num)
  burst.time<-c(start.times=bursts[min(clust), "start"], end.times=bursts[max(clust), "end"] )
  }
  bursts.comb<-sapply(unique(bursts$clusid), combine.clusters)
  start.times<-bursts.comb[1,]
  end.times<-bursts.comb[2,]
  N.burst<-length(start.times)
  beg<-sapply(start.times, function(x) which.min(abs(spike.train-x)))
  end<- sapply(end.times, function(x) which.min(abs(spike.train-x)))
  IBI<-c(NA, start.times[-1]-end.times[-N.burst])
  len<-end-beg+1
  durn<-end.times-start.times
  mean.isis<-durn/(len-1)
  result<-cbind(beg=beg, end=end, IBI=IBI, len=len, durn=durn, mean.isis=mean.isis, SI=rep(thresh, N.burst))
  result
}