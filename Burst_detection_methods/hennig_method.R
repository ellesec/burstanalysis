hennig.method<-function(st, cutoff.prob=0.05) {
  bursts<-NULL
  result<-NULL
  allisi<-diff(st)
  if (length(allisi)<1) {
    result<- NA
  } else {
  isi.rank<-rank(allisi) #1 is smallest ISI
  st.length<-ceiling(max(st))
  spike.counts<-NULL
  for (i in 0:st.length-1) {
    spike.counts[i+1]<-sum((st>=i)*(st<(i+1)))
  }
  sc.hist<-hist(spike.counts, nclass=200)
  p.dist<-1-cumsum(sc.hist$counts/sum(sc.hist$counts))
  cutoff.indx<-sum(p.dist>cutoff.prob)
  theta.c<-max(c(2, ceiling(sc.hist$mids[cutoff.indx])))
  theta.c.off<-theta.c*0.5
  isi.rel.rank<-isi.rank/max(isi.rank)
  
  #t<-st[2:length(st)]#Why two??
  t<-st
  j<-1
  burst.on<-0
  bc<-1
  dt<-1
  burst.time<-NULL
  burst.end<-NULL
  burst.dur<-NULL
  burst.size<-NULL
  burst.beg<-NULL
  while (j<length(allisi)-theta.c) {
    if (burst.on==0 && isi.rel.rank[j]<0.5) {
      if (t[j+theta.c]<t[j]+dt){
        burst.on<-1
        burst.time[bc]<-t[j]
        burst.beg[bc]<-j
        brc<-j
      }
    } else if (burst.on==1) { 
      if (t[j+theta.c.off]>t[j]+dt) {
        burst.end[bc]<-t[j]
        burst.dur[bc]<-t[j]-burst.time[bc]
        burst.size[bc]<-j-brc
        bc<-bc+1
        burst.on<-0
      }
    }
    j<-j+1
  }
  if (burst.on==1) {
    tmp<-t[j]-burst.time[bc]
    burst.end[bc]<-burst.time[bc]+tmp
    burst.dur[bc]<-t[j]-burst.time[bc]
    burst.size[bc]<-j-brc
    bc<-bc+1
  }
  N.burst<-length(burst.time)
  if (N.burst<1) {
    result<-NA
  } else {
  end<-burst.beg+burst.size
  IBI<-c(NA, burst.time[-1]-burst.end[-N.burst])
  len<-burst.size+1
  mean.isis<-burst.dur/(len-1)
  result<-cbind(beg=burst.beg, end=end, IBI=IBI, len=len, durn=burst.dur, mean.isis=mean.isis, SI=rep(1, N.burst))
  }
  }
  result
}