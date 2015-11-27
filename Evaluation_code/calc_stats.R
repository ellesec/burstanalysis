#Calculate burst statistics from running method on a spike train. Input is spike train, method type, and threshold levels
#Output is statistics of bursting behaviours
run.method<-function(spike.train, method, cutoff=NULL) {
  bursts <- switch(method, mi = MI.method(spike.train), 
                  ps = PS.method(spike.train, cutoff),
                  rs = RS.method(spike.train, cutoff), rgs=RGS.method(spike.train), hsmm=HSMM.method(spike.train, cutoff), 
                  cma=CMA.method(spike.train, brs.incl=cutoff), hennig=hennig.method(spike.train),  logisi.pasq=logisi.pasq.method(spike.train, cutoff), stop(method, " : no such method for burst analysis"))
if (class(bursts)=="list") {
 burst.stats<- lapply(bursts, function(x) calc.st.stats(x, spike.train))
 names(burst.stats)<-cutoff
} else {
  burst.stats<-calc.st.stats(bursts, spike.train)
}
  
burst.stats
}

#Calculates statistics based on detected bursts and spike train
calc.st.stats<-function(bursts, spike.train, duration=300) {
  channels<-1
  spikes<-length(spike.train)
  mean.freq<-round(spikes/duration,3)
  nbursts<-max(dim(bursts)[1],0)
  bursts.per.sec<-max(round(nbursts/duration,3),0)
  bursts.per.min<-bursts.per.sec*60
  durations<-get.details(bursts, "durn")
  mean.dur<-mean(durations)
  sd.dur<-sd(durations)
  ISIs<-diff(spike.train)
  mean.ISIs<-mean(ISIs)
  ns<-get.details(bursts, "len")
  mean.spikes<-round(mean(ns),3)
  sd.spikes<-round(sd(ns),3)
  total.spikes.in.burst<-sum(ns)
  per.spikes.in.burst<-round(100*(total.spikes.in.burst/spikes), 3)
  IBIs<-get.details(bursts, "IBI")
  mean.IBIs<-mean(IBIs, na.rm=TRUE)
  sd.IBIs<-sd(IBIs, na.rm=TRUE)
  cv.IBIs<-round(sd.IBIs/mean.IBIs, 3)
  cv.IBIs<-ifelse(is.infinite(cv.IBIs), NA, cv.IBIs)
  df <- data.frame( spikes = spikes, mean.freq = mean.freq, 
                   nbursts = nbursts, bursts.per.sec = bursts.per.sec, bursts.per.min = bursts.per.min, 
                    mean.dur = round(mean.dur,3),  
                   mean.spikes = mean.spikes, per.spikes.in.burst = per.spikes.in.burst, 
                   per.spikes.out.burst = round(100 - per.spikes.in.burst, 
                                                3),  mean.isis = mean.ISIs, 
                   mean.IBIs = mean.IBIs,  cv.IBIs = cv.IBIs)
}

#Gets values from burst matrix
get.details<-function(b, index) {
if (length(b) > 1) {
  b[, index]
}
else {
  0
}
}