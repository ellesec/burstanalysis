#Compares detected bursts with 'ground truth' behaviour. Input is spike train, method type, true beginning and true end of bursts and threshold value
#Output is statistics of bursting behaviours
eval.method<-function(spike.train, method, beg.true, end.true, cutoff=NULL) {
  bursts <- switch(method, mi = MI.method(spike.train), 
                   ps = PS.method(spike.train, cutoff),
                   rs = RS.method(spike.train, cutoff), rgs=RGS.method(spike.train, thresh=cutoff), hsmm=HSMM.method(spike.train, cutoff),
                   cma=CMA.method(spike.train), hennig=hennig.method(spike.train, cutoff), logisi.pasq=logisi.pasq.method(spike.train, cutoff), stop(method, " : no such method for burst analysis"))
  if (class(bursts)=="list") {
    result<-lapply(bursts, function(x) eval.stats(spike.train, beg.true, end.true, x))
  } else {
    result<-eval.stats(spike.train, beg.true, end.true, bursts)
  }
  
result
}


eval.stats<-function(spike.train, beg.true, end.true, bursts) {
  n.spikes<-length(spike.train)
  if (is.null(beg.true)) {
    #If no true bursts, and method finds no bursts, true positive rate =1, false positive rate NA
    if (is.na(bursts[1])) {
      return(c(1,NA))
    }
    false.pos.rate<-sum(bursts[,"len"])/n.spikes
    #If no true bursts, and method finds bursts, true positive rate NA, false positive rate is number of spikes found in bursts
    res<-c(NA, false.pos.rate)
  } else{
    #If method finds no bursts when true bursts do exist, both rates are zero
    if (is.na(bursts[1])) {
      return(c(0,0))
    }
  n.bursts<-dim(bursts)[1]
  method.bursts<-rep(0, n.spikes)
 method.indx<-unlist(sapply(1:n.bursts, function(x) seq(bursts[x,"beg"],bursts[x,"end"])))
 method.bursts[method.indx]<-1
 true.bursts<-rep(0, n.spikes)
 true.indx<-unlist(sapply(1:length(beg.true), function(x) seq(beg.true[x], end.true[x])))
 true.bursts[true.indx]<-2
 true.burst.length<-length(true.indx)
 diff.bursts<-true.bursts-method.bursts
 true.pos<-which(diff.bursts>(1-.Machine$double.eps) & diff.bursts<(1+.Machine$double.eps))
 false.pos<-which(diff.bursts< (-1*.Machine$double.eps))
 res<-c(length(true.pos)/true.burst.length, length(false.pos)/(n.spikes-true.burst.length))
  }
}


evaluate <- function(spike.times, start.true, end.true, start.method, end.method){
  

  
    # Use spike.times and true start/end times to label spikes burst (=1) vs. no burst (=0)
    true.ID <- rep(0, length(spike.times))
    for (x in 1:length(start.true)){
      burst <- which(spike.times <= end.true[x] & spike.times >= start.true[x])
      true.ID[burst] <- 1
    }
    
    # Use spike.times and method start/end times to label spikes burst vs. no burst
    method.ID <- rep(0, length(spike.times))
    for (x in 1:length(start.method)){
      burst <- which(spike.times <= end.method[x] & spike.times >= start.method[x])
      method.ID[burst] <- 1
    }
    
    # Compare method.ID against true.ID
    total.pos <- length(which(true.ID==1))
    total.neg <- length(which(true.ID==0))
    score <-(2*true.ID)-method.ID       #true.pos=1, true.neg=0, false.pos=-1, false.neg=2
    true.pos <- length(which(score==1))
    false.pos <- length(which(score==-1))
    true.pos.rate <- true.pos/total.pos
    false.pos.rate <- false.pos/total.neg
    no.bursts <- length(start.method)
  return(list(true.pos.rate, false.pos.rate, no.bursts))
}
