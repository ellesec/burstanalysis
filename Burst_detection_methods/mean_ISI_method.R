MISI.method<-function(spike.train) {
isis<-diff(spike.train)
mn.isi<-mean(isis)
L<-isis[isis<mn.isi]
ML<-mean(L)
N<-length(isis)
bursts<-NULL
i<-1
while(i<N){
  cum.avg<-cumsum(isis[i:N]) / seq_along(isis[i:N])
  joint.isis<-(cum.avg<=ML)[-1]
  if (any(joint.isis)){
  end.b<-max(which(joint.isis))+i
  bursts<-rbind(bursts, c(i, end.b))
  i<-end.b+1
  } else{
    i<-i+1
  }
}
N.burst<-dim(bursts)[1]
beg<-bursts[,1]
end<-bursts[,2]+1
start.times<-spike.train[beg]
end.times<-spike.train[end]
IBI<-c(NA, start.times[-1]-end.times[-N.burst])
len<-end-beg+1
durn<-end.times-start.times
mean.isis<-durn/(len-1)
result<-cbind(beg=beg, end=end, IBI=IBI, len=len, durn=durn, mean.isis=mean.isis, SI=rep(1, N.burst))
result
}