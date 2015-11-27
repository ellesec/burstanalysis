library(sjemea)
library(e1071)

##Function to calculate and plot bursts, based on CMA method (Kapucu et al., 2012).
##Input is spike train. If brs.incl is set to true, method will include burst related
##spikes. min.val is the minumum number of spikes on a spikes train for the method to be run.
##If plot = TRUE, spike train with bursts labelled is plotted.
CMA.method<-function(spike.train, brs.incl=TRUE, min.val=3, plot=FALSE) {
  #Do not perform burst detection if less than min.val spikes in the spike train
  if (length(spike.train)<min.val) {
    result<-NA
    return(result)
  }
  isi<-diff(spike.train)
  isi.range<-max(isi)-min(isi)
  eps<-isi.range/1000
  if (isi.range<0.001){
    breaks1<-seq(0, max(isi)+isi.range/10, isi.range/10)
    hist.isi<-hist(isi, breaks =breaks1, plot=FALSE) #if ISI range very small (<0.001), use smaller bins
  } else {
    hist.isi<-hist(isi, breaks =seq(0, max(isi)+eps, eps), plot=FALSE) #Create histogram with approx 1000 bins
  }
  CMA<-cumsum(hist.isi$counts)/seq(1,length(hist.isi$counts),1)
  CMAm<-max(CMA)
  m<-min(which(CMA==CMAm))
  alpha.values<-data.frame(max=c(1, 4, 9, 1000), a1=c(1, 0.7, 0.5, 0.3), a2=c(0.5, 0.5, 0.3,0.1)) #Alpha value scale
  skew<-skewness(CMA)
  if (is.na(skew)){
    result<-NA
    return(result)
  }
  diff.skew<-alpha.values[,"max"]-skew
  alpha.indx<-which(diff.skew==min(diff.skew[diff.skew>0]))
  alpha1<-alpha.values[alpha.indx, "a1"]
  alpha2<-alpha.values[alpha.indx, "a2"]
  cutoff<-which.min(abs(CMA[m:length(CMA)] - (alpha1*CMAm)))+(m-1) #Cutoff set at bin closest in value to alpha1*CMAm
  xt<-hist.isi$mids[cutoff] #maxISI
  cutoff2<-which.min(abs(CMA[m:length(CMA)] - (alpha2*CMAm)))+(m-1) #Burst related spikes cutoff set at bin closest in value to alpha2*CMAm
  xt2<-hist.isi$mids[cutoff2] #maxISI for burst related spikes
  bursts<-find.bursts(spike.train,xt) #Find burst cores
  #If brs.incl=TRUE, then extend bursts to include burst related spikes
  if (brs.incl && !is.null(dim(bursts)[1])){
    brs<-find.bursts(spike.train, xt2)
    burst.adj<-NULL
    for (i in 1:dim(bursts)[1]){
    burst.between<-apply(brs[,1:2], 1, function(x) between.bursts(bursts[i,1:2], x)) #Find burst related spikes which surround bursts
    which.between<-which(unlist(sapply(burst.between, function(x) sum(!is.na(x))))>0)
    if (which.between) {
      burst.adj<-rbind(burst.adj, burst.between[[which.between]])
    } else {
      burst.adj<-rbind(burst.adj, bursts[i, 1:2])
    }
  }
  
  burst.adj<-unique(burst.adj) #Remove any repeated bursts
  N<-dim(burst.adj)[1]
  beg<-burst.adj[,1]
  end<-burst.adj[,2]
  ibi<-c(NA, spike.train[beg[-1]]-spike.train[end[-N]])
  len<-end-beg+1
  durn<-spike.train[end]-spike.train[beg]
  bursts<-cbind(beg=beg, end=end, IBI=ibi, len=len, durn=durn, mean.isis=durn/len, SI=1)
  }

  if (is.null(dim(bursts)[1])){
    bursts<-NA
  }
  bursts
}
  
    
    between.bursts<-function(burst1, burst2) {
      if (burst2[1]<=burst1[1] & burst2[2]>=burst1[2]) {
        burst<-c(min(burst1[1], burst2[1]), max(burst2[1], burst2[2]))
      } else {
        burst<-NA
      }
      burst
    }
    
    #Add back in IBIs

find.bursts<-function(spike.train, xt){
isi<-diff(spike.train)
indxs<-which(isi<xt)
burst.breaks<- c(0, which(diff(indxs) >1), length(indxs))
isi.list<-sapply(seq(length(burst.breaks) - 1), function(i) indxs[(burst.breaks[i] + 1):burst.breaks[i+1]])
burst.indx<-which(sapply(isi.list, length)>1)
if(length(burst.indx)){
beg<-sapply(isi.list[burst.indx], function(x) min(x))
end<-sapply(isi.list[burst.indx], function(x) max(x))+1
N<-length(beg)
ibi<-c(NA, spike.train[beg[-1]]-spike.train[end[-N]])
len<-end-beg+1
durn<-spike.train[end]-spike.train[beg]
bursts<-cbind(beg=beg, end=end, IBI=ibi, len=len, durn=durn, mean.isis=durn/(len-1), SI=1)
} else {
  bursts<-NA
}
bursts
}

