library(sjemea)
library(e1071)

##Function to calculate and plot bursts, based on CMA method (Kapucu et al., 2012).
##Input is spike train. If brs.incl is set to true, method will include burst related
##spikes. min.val is the minumum number of spikes on a spikes train for the method to be run.
##If plot = TRUE, spike train with bursts labelled is plotted.
CMA.method<-function(s1, brs.incl=TRUE, min.val=3, plot=FALSE) {
  if (length(s1)<min.val) {
    result<-NA
    return(result)
  }
  s1.isi<-diff(s1)
  isi.range<-max(s1.isi)-min(s1.isi)
  eps<-isi.range/1000
  if (isi.range<0.001){
    breaks1<-seq(0, max(s1.isi)+isi.range/10, isi.range/10)
    hist.isi<-hist(s1.isi, breaks =breaks1, plot=FALSE) #if ISI range very small (<0.001), use smaller bins
  } else {
    hist.isi<-hist(s1.isi, breaks =seq(0, max(s1.isi)+eps, isi.range/1000), plot=FALSE) #
  }
  CMA<-cumsum(hist.isi$counts)/seq(1,length(hist.isi$counts),1)
  CMAm<-max(CMA)
  m<-min(which(CMA==CMAm))
  alpha.values<-data.frame(max=c(1, 4, 9, 1000), a1=c(1, 0.7, 0.5, 0.3), a2=c(0.5, 0.5, 0.3,0.1))
  skew<-skewness(CMA)
  if (is.na(skew)){
    result<-NA
    return(result)
  }
  diff.skew<-alpha.values[,"max"]-skew
  alpha.indx<-which(diff.skew==min(diff.skew[diff.skew>0]))
  alpha1<-alpha.values[alpha.indx, "a1"]
  alpha2<-alpha.values[alpha.indx, "a2"]
  cutoff<-which.min(abs(CMA[m:length(CMA)] - (alpha1*CMAm)))+(m-1) 
  xt<-hist.isi$mids[cutoff]
  cutoff2<-which.min(abs(CMA[m:length(CMA)] - (alpha2*CMAm)))+(m-1) 
  xt2<-hist.isi$mids[cutoff2]
  st<-cumsum(s1.isi)
  bursts<-CMA.find.bursts(s1,xt)
  brs<-CMA.find.bursts(s1, xt2)
  
  if(dim(bursts)[1]) {
    burst.adj<-data.frame(beg=rep(0, dim(bursts)[1]), end=rep(0, dim(bursts)[1]) )
    
    if (brs.incl & dim(brs)[1]) {
      for (i in 1:dim(bursts)[1]) {
        for (j in 1:dim(brs)[1]) {
          if(is.between(bursts[i,1], brs[j,1], brs[j,2]) | is.between(bursts[i,2], brs[j,1], brs[j,2]))
          {
            burst.adj$beg[i]<-min(bursts[i,1], brs[j,1])
            burst.adj$end[i]<-max(bursts[i,2], brs[j,2])
            break
          } else {
            burst.adj$beg[i]<-bursts[i,1]
            burst.adj$end[i]<-bursts[i,2]
          }
          if(brs[j,2]>bursts[i,2]) {
            break
          }
        }
      }
      
      
      diff.begs<-diff(burst.adj[,"beg"])
      rep.bursts.begs<-which(diff.begs==0)
      if (any(rep.bursts.begs)) {
        burst.adj<-burst.adj[-rep.bursts.begs,]
      }
      diff.ends<-diff(burst.adj[,"end"])
      rep.bursts.end<-which(diff.ends==0)+1
      if (any(rep.bursts.end)) {
        burst.adj<-burst.adj[-rep.bursts.end,]
      }
      
    } else {
      burst.adj<-data.frame(beg = bursts[,1], end = bursts[,2])
    }
    
    n<-2
    IBI<-c(NA)
    while (n<=dim(burst.adj)[1]) {
      IBI<-c(IBI, s1[burst.adj[n,1]]-s1[burst.adj[n-1,2]])
      n<-n+1
    }
    burst.adj<-cbind(burst.adj, IBI)
    names(burst.adj)<-c("beg", "end", "IBI")
    
    
    is.between<-function(x,a,b){ betw<-0
                                 if(x>=a &x<=b)
                                 {
                                   betw<-1
                                 }
                                 betw
    }
    
    #Add back in IBIs

    merge.bursts <- which(burst.adj[,"IBI"] < xt2)
    if (any(merge.bursts)) {
      for (n in rev(merge.bursts)) {
        burst.adj[n-1, "end"]<- burst.adj[n, "end"]
      }
      burst.adk <- burst.adj[-merge.bursts,,drop=FALSE] 
      burst.adj<-burst.adk
    }
   durn<-s1[burst.adj$end]-s1[burst.adj$beg]
   len<-burst.adj$end-burst.adj$beg+1
   mean.isis<-durn/(len-1)
   result<-cbind(beg=burst.adj$beg, end=burst.adj$end, IBI=burst.adj$IBI, len=len, durn=durn, mean.isis=mean.isis, SI=1)
    
  } else {
    result<-NA
  }
  result
}


#Find bursts using defined cutoff 
CMA.find.bursts <- function(spikes, xt,debug=FALSE) {
  
  ## For one spike train, find the burst using max interval method.
  ## e.g.
  
  nspikes = length(spikes)
  
  no.bursts = matrix(nrow=0,ncol=1)  
  
  ## Create a temp array for the storage of the bursts.  Assume that
  ## it will not be longer than Nspikes/2 since we need at least two
  ## spikes to be in a burst.
  
  max.bursts <- floor(nspikes/2)
  bursts <- matrix(NA, nrow=max.bursts, ncol=3)
  colnames(bursts) = c("beg", "end", "IBI")
  burst <- 0                            #current burst number
  
  ## Phase 1 -- burst detection.  Here a burst is defined as starting
  ## when two consecutive spikes have an ISI less than BEG.ISI apart.
  ## The end of the burst is given when two spikes have an ISI greater
  ## than END.ISI.
  
  ## Find ISIs closer than beg.isi, and end with end.isi.
  
  
  ## LAST.END is the time of the last spike in the previous burst.
  ## This is used to calculate the IBI.
  ## For the first burst, this is no previous IBI
  last.end = NA;                        #for first burst, there is no IBI.
  
  n = 2
  in.burst = FALSE
  
  
  while ( n <= nspikes) {
    
    next.isi = spikes[n] - spikes[n-1]
    if (in.burst) {
      if (next.isi > xt) {
        ## end of burst
        end = n-1; in.burst = FALSE
        
        
        ibi =  spikes[beg] - last.end; last.end = spikes[end]
        res = c(beg, end, ibi)
        burst = burst + 1
        if (burst > max.bursts) {
          print("too many bursts!!!")
          browser()
        }
        bursts[burst,] <- res
      }
    } else {
      ## not yet in burst.
      if (next.isi <= xt) {
        ## Found the start of a new burst.
        beg = n-1; in.burst = TRUE
      }
    }
    n = n+1
  }
  
  ## At the end of the burst, check if we were in a burst when the
  ## train finished.
  if (in.burst) {
    end = nspikes
    ibi =  spikes[beg] - last.end
    res = c(beg, end, ibi)
    burst = burst + 1
    if (burst > max.bursts) {
      print("too many bursts!!!")
      browser()
    }
    bursts[burst,] <- res
  }
  
  ## Check if any bursts were found.
  if (burst > 0 ) {
    ## truncate to right length, as bursts will typically be very long.
    bursts = bursts[1:burst,,drop=FALSE]
  } else {
    ## no bursts were found, so return an empty structure.
    return(no.bursts)
  }
  
  if (debug) {
    print("End of phase1\n")
    print(bursts)
  }
  
  
  ## Phase 2 -- merging of bursts.  Here we see if any pair of bursts
  ## have an IBI less than MIN.IBI; if so, we then merge the bursts.
  ## We specifically need to check when say three bursts are merged
  ## into one.
  
  
  ## Phase 3 -- remove small bursts: less than min duration (MIN.DURN), or
  ## having too few spikes (less than MIN.SPIKES).
  ## In this phase we have the possibility of deleting all spikes.
  
  ## LEN = number of spikes in a burst.
  ## DURN = duration of burst.
  len = bursts[,"end"] - bursts[,"beg"] + 1
  durn = spikes[bursts[,"end"]] - spikes[bursts[,"beg"]]
  bursts = cbind(bursts, len, durn)
  
  rejects = which ( ( len < 3) )
  
  if (any(rejects)) {
    bursts = bursts[-rejects,,drop=FALSE]
  }
  
  if (nrow(bursts) == 0) {
    ## All the bursts were removed during phase 3.
    bursts = no.bursts
  } else {
    ## Compute mean ISIS
    len = bursts[,"end"] - bursts[,"beg"] + 1
    durn = spikes[bursts[,"end"]] - spikes[bursts[,"beg"]]
    mean.isis = durn/(len-1)
    
    ## Recompute IBI (only needed if phase 3 deleted some cells).
    if (nrow(bursts)>1) {
      ibi2 = c(NA, calc.ibi(spikes, bursts))
    } else {
      ibi2 = NA
    }
    bursts[,"IBI"] = ibi2
    
    SI = rep(1, length(mean.isis ))
    bursts = cbind(bursts, mean.isis, SI)
  }
  
  ## End -- return burst structure.
  bursts
  
}


ibis = burst.adj[,"IBI"]
merge.bursts = which(ibis < xt2)

if (any(merge.bursts)) {
  ## Merge bursts efficiently.  Work backwards through the list, and
  ## then delete the merged lines afterwards.  This works when we
  ## have say 3+ consecutive bursts that merge into one.
  
  for (burst in rev(merge.bursts)) {
    burst.adj[burst-1, "end"] = burst.adj[burst, "end"]
    burst.adj[burst, "end"] = NA         #not needed, but helpful.
  }
  burst.adk = burst.adj[-merge.bursts,,drop=FALSE] #delete the unwanted info.
}

if (debug) {
  print("End of phase 2\n")
  print(bursts)
}

