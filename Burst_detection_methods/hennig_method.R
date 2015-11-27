#ISI Rank threshold method
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
    spike.counts[i+1]<-sum((st>=i)*(st<(i+1))) #calculate spike count of 1s intervals
  }
  sc.hist<-hist(spike.counts, nclass=200, plot=FALSE)
  p.dist<-1-cumsum(sc.hist$counts/sum(sc.hist$counts)) 
  cutoff.indx<-sum(p.dist>cutoff.prob)
  theta.c<-max(c(2, ceiling(sc.hist$mids[cutoff.indx]))) #set theta_C to value where probability of spike counts is equal to cutoff.index (default 0.05)
  theta.c.end<-theta.c*0.5 #cutoff to end a burst
  isi.rel.rank<-isi.rank/max(isi.rank) #calculate relative rank of each isi
  
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
    if (burst.on==0 && isi.rel.rank[j]<0.5) { #burst begins when rank of isi<0.5
      if (t[j+theta.c]<t[j]+dt){ 
        burst.on<-1
        burst.time[bc]<-t[j]
        burst.beg[bc]<-j
        brc<-j
      }
    } else if (burst.on==1) { 
      if (t[j+theta.c.end]>t[j]+dt) {
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