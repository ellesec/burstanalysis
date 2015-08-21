#Long bursts
pois.burst.long<-function(rate, period, start) {
  burst.times<-sort(runif(rpois(1,period*rate), min=start, max=start+period))
  no.bursts<-length(burst.times)
  burst.length<-rpois(no.bursts, mean.length)
  burst.pos<-sapply(burst.length, function(x) runif(x, min=-1.5, max=1.5))
  burst.pos.sorted<-lapply( burst.pos, sort)
  spike.times<-list()
  for ( i in 1:length(burst.times)){
    spike.times[[i]]<-burst.pos.sorted[[i]]+burst.times[i]
  }
  spike.times2<-spike.times[which(burst.length>2)]
  burst.beg<-sapply(spike.times2, min)
  burst.end<-sapply(spike.times2, max)
  N<-length(burst.beg)
  overlap<-which((burst.beg[-1]-burst.end[-N])<0.5)
  if(any(overlap)){
    spike.times2<-spike.times2[-overlap]
    burst.beg<-burst.beg[-overlap]
    burst.end<-burst.end[-overlap]
  }
  sts<-sort(unlist(spike.times2))
  sts<-sts[which(sts>0)]
  N.bursts<-length(burst.beg)
  insim<-list()
  insim$spks<-sts
  insim$burst.beg<-burst.beg
  insim$burst.end<-burst.end
  insim$num.bursts<-length(burst.beg)
  insim
}



#High frequency
pois.burst.high<-function(rate, period, start) {
  burst.times<-sort(runif(rpois(1,period*rate), min=start, max=start+period))
  no.bursts<-length(burst.times)
  burst.length<-rpois(no.bursts, mean.length)
  burst.pos<-sapply(burst.length, function(x) runif(x, min=-0.25, max=0.25))
  burst.pos.sorted<-lapply( burst.pos, sort)
  spike.times<-list()
  for ( i in 1:length(burst.times)){
    spike.times[[i]]<-burst.pos.sorted[[i]]+burst.times[i]
  }
  spike.times2<-spike.times[which(burst.length>2)]
  burst.beg<-sapply(spike.times2, min)
  burst.end<-sapply(spike.times2, max)
  sts<-sort(unlist(spike.times2))
  sts<-sts[which(sts>0)]
  N.bursts<-length(burst.beg)
  insim<-list()
  insim$spks<-sts
  insim$burst.beg<-burst.beg
  insim$burst.end<-burst.end
  insim$num.bursts<-length(burst.beg)
  insim
}

#Bursts with noise
pois.burst.noisy<-function(rate, period, start) {
  burst.times<-sort(runif(rpois(1,period*rate), min=start, max=start+period))
  no.bursts<-length(burst.times)
  #burst.length<-NULL
  #for (i in 1:no.bursts) {
  #  mn<- round(burst.times[i], -2)/100+4
  #  burst.length[i]<-rpois(1, mn)
  #}
  burst.length<-rpois(no.bursts, mean.length)
  burst.pos<-sapply(burst.length, function(x) runif(x, min=-0.4, max=0.4))
  burst.pos.sorted<-lapply( burst.pos, sort)
  spike.times<-list()
  for ( i in 1:length(burst.times)){
    spike.times[[i]]<-burst.pos.sorted[[i]]+burst.times[i]
  }
  spike.times2<-spike.times[which(burst.length>2)]
  while (spike.times2[[1]][1]<0){
    spike.times2<-spike.times2[-1]
  }
  burst.beg<-sapply(spike.times2, min)
  burst.end<-sapply(spike.times2, max)
  overlap<-burst.beg[-1]-burst.end[-length(burst.end)]
  ol<-which(overlap<0.5)
  if (length(ol)>0) {
    spike.times2<-spike.times2[-ol]
    burst.beg<-burst.beg[-ol]
    burst.end<-burst.end[-ol]
  }
  sts<-sort(unlist(spike.times2))
  sts<-sts[which(sts>0)]
  
  
  
  
  gst<-rgamma(200, 1,  rate=0.5)
  q<-quantile(gst, 0.1)
  gst<-gst[-which(gst<q)]
  gsum<-cumsum(gst)
  noise.filtered<-gsum[which(gsum<period)]
  
  burst.mids<-(burst.beg-burst.end)/2+burst.beg
  
  
  
  spk.rem<-NULL
  for (i in 1:length(burst.mids)) {
    burst.isis<-abs(noise.filtered-burst.mids[i])
    spk.rem<-c(spk.rem, which(burst.isis<0.9))
  }
  noise.filtered2<-noise.filtered[-spk.rem]
  
  spks<-sort(c(noise.filtered2, sts))
  insim<-list()
  insim$spks<-spks
  insim$burst.beg<-sapply(burst.beg, function(x) which(spks==x))
  insim$burst.end<-sapply(burst.end, function(x) which(spks==x))
  insim
}

#Regular bursts
pois.burst.noiseless<-function(rate, period, start) {
  burst.times<-sort(runif(rpois(1,period*rate), min=start, max=start+period))
  no.bursts<-length(burst.times)
  burst.length<-rpois(no.bursts, mean.length)
  burst.pos<-sapply(burst.length, function(x) runif(x, min=-0.15, max=0.15))
  burst.pos.sorted<-lapply( burst.pos, sort)
  spike.times<-list()
  for ( i in 1:length(burst.times)){
    spike.times[[i]]<-burst.pos.sorted[[i]]+burst.times[i]
  }
  spike.times2<-spike.times[which(burst.length>2)]
  burst.beg<-sapply(spike.times2, min)
  burst.end<-sapply(spike.times2, max)
  overlap<-burst.beg[-1]-burst.end[-length(burst.end)]
  ol<-which(overlap<0.5)
  if (length(ol)>0) {
    spike.times2<-spike.times2[-ol]
  }
  burst.beg<-sapply(spike.times2, min)
  burst.end<-sapply(spike.times2, max)
  sts<-sort(unlist(spike.times2))
  sts<-sts[which(sts>0)]
  
  insim<-list()
  insim$spks<-sts
  insim$burst.beg<-burst.beg
  insim$burst.end<-burst.end
  insim$num.bursts<-length(burst.beg)
  insim
}