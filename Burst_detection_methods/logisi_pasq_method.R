library(pracma)

#Function to run logISI method
logisi.pasq.method<-function(spike.train, cutoff=0.1){
  cutoff<-ifelse(is.null(cutoff), 0.1, cutoff)
  if (length(spike.train)>3) {
      isi.low <- logisi.break.calc(spike.train, cutoff) #Calculates threshold as isi.low
    if (is.null(isi.low) || isi.low>=1 ){
      logisi.par <- list(min.ibi=0,   min.durn=0, min.spikes=3,
      isi.low=cutoff) #If no value for isi.low found, or isi.low above 1 second, find bursts using threshold equal to cutoff (default 100ms)
      result<-logisi.find.burst(spike.train, logisi.par)
    } else if (isi.low<0) {
      result<-NA
    } else if (isi.low>cutoff & isi.low <1) {
      logisi.par <- list(min.ibi=isi.low,   min.durn=0, min.spikes=3,
                         isi.low=cutoff) #If isi.low >cutoff, find bursts using threshold equal to cutoff (default 100ms)
      bursts<-logisi.find.burst(spike.train, logisi.par)
      if (!is.na(bursts)[1]){
        logisi.par2 <- list(min.ibi=0,   min.durn=0, min.spikes=3,
        isi.low=isi.low) #If bursts have been found, add burst related spikes using threshold of isi.low
        brs<-logisi.find.burst(spike.train, logisi.par2)
        result<-add.brs(bursts, brs, spike.train)
      } else {
        result<-bursts
      }
    } else {
      logisi.par <- list(min.ibi=0,   min.durn=0, min.spikes=3,
      isi.low=isi.low) #If isi.low<cutoff, find bursts using a threshold equal to isi.low
      result<-logisi.find.burst(spike.train, logisi.par)
    }
    
  } else {
    result<-NA
  }
  result
}

#Finds peaks in logISI histogram
get.peaks<-function(h, Pd=2, Th=0, Np=NULL){
  m<-0
  L<-length(h$density)
  j<-0
  Np<-ifelse(is.null(Np), L, Np)
  pks<-NULL
  locs<-NULL
  void.th<-0.7
  while((j<L)&&(m<Np)){
    j<-j+1
    endL<-max(1,j-Pd)
    if (m>0 && j<min(c(locs[m]+Pd, L-1))){
      j<-min(c(locs[m]+Pd, L-1))
      endL<-j-Pd
    }
    endR<-min(L, j+Pd)
    temp<-h$density[endL:endR]
    aa<-which(j==endL:endR)
    temp[aa]<--Inf
    if (Pd>1){
      idx1<-max(1, aa-2)
      idx2<-min(aa+2, length(temp))
      idx3<-max(1, aa-1)
      idx4<-min(aa+1, length(temp))
      if (sum((h$density[j]>(temp[c(1:idx1, idx2:length(temp))]+Th))==FALSE)==0 && sum((h$density[j]>(temp[idx3:idx4]))==FALSE)==0 && j!=1 && j!=L){
        m<-m+1
        pks[m]<-h$density[j]
        locs[m]<-j
      } } else if (sum((h$density[j]>(temp+Th))==FALSE)==0 ) {
        m<-m+1
        pks[m]<-h$density[j]
        locs[m]<-j
      }
    
  }
  ret<-data.frame(pks=pks, locs=locs)
}

#Function to find cutoff threshold.
find.thresh<-function(h, ISITh=100){
  void.th<-0.7
  gp<-get.peaks(h)
  num.peaks<-length(gp$pks)
  pkx<-h$breaks[gp$locs]
  intra.indx<-which(pkx<ISITh)
  if(length(intra.indx)>=1){
    max.intra<-max(gp$pks[intra.indx])
    max.idx<-which.max(gp$pks[intra.indx])
  } else {
    return(-1000)
  }
  x1<-pkx[max.idx]
  y1<-max.intra
  locs1<-gp$locs[max.idx]
  num.peaks.after.burst<-num.peaks-max.idx
  if (num.peaks.after.burst==0){
    return(NULL)
  } else {
    gp2<-gp[(max.idx+1):num.peaks,]
    ymin<-sapply(gp2$locs, function(x) min(h$density[locs1:x]))
    xmin<-sapply(gp2$locs, function(x) which.min(h$density[locs1:x]))+locs1-1
    voidParameter<-1-(ymin/sqrt(y1*gp2$pks))
  }
  indxvoid<-suppressWarnings(min(which(voidParameter>=void.th)))
  if(is.infinite(indxvoid)) {
    flags<-c(1,0)
    return(NULL)
  } else {
    ISImax<-h$breaks[xmin[indxvoid]]
    return(ISImax)
  }
}

#Calculates cutoff for burst detection
logisi.break.calc<-function(st, cutoff){
  isi<-diff(st)*1000
  max.isi<-ceiling(log10(max(isi)))
  isi<-isi[isi>=1]
  br<-logspace(0, max.isi, 10*max.isi)
  h<-hist(isi, breaks=br, plot=FALSE)
  h$density<-h$counts/sum(h$counts)
  h$density<-lowess(h$density, f=0.05)$y
  thr<-find.thresh(h, cutoff*1000)
  if(!is.null(thr)){
    thr<-thr/1000
  }
  thr
}




###Function to add burst related spikes to edges of bursts
add.brs<-function(bursts, brs, spike.train){
is.between<-function(x,a,b){ betw<-0
                             if(x>=a &x<=b)
                             {
                               betw<-1
                             }
                             betw
}
                             
burst.adj<-data.frame(beg=rep(0, dim(bursts)[1]), end=rep(0, dim(bursts)[1]) )
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
start.times<-spike.train[burst.adj$beg]
end.times<- spike.train[burst.adj$end]
durn<-end.times-start.times
len<-burst.adj$end-burst.adj$beg+1
mean.isis<-durn/(len-1)
N.burst<-dim(burst.adj)[1]
IBI<-c(NA, start.times[-1]-end.times[-N.burst])
result<-cbind(beg=burst.adj$beg, end=burst.adj$end, IBI=IBI, len=len, durn=durn, mean.isis=mean.isis, SI=rep(1, N.burst))
result

}


##Function for finding bursts, taken from sjemea
logisi.find.burst<- function(spikes, par, debug=FALSE) {
  
  ## For one spike train, find the burst using log isi method.
  ## e.g.
  ## find.bursts(s$spikes[[5]])
  ## init.
  ## params currently in LOGISI.PAR
  ##
  
  no.bursts = NA;                       #value to return if no bursts found.
  
  
  ##beg.isi =    par$beg.isi
  ##end.isi =    par$end.isi
  min.ibi =      par$min.ibi
  min.durn =     par$min.durn
  min.spikes =   par$min.spikes
  isi.low =      par$isi.low
  
  nspikes = length(spikes)
  
  ## Create a temp array for the storage of the bursts.  Assume that
  ## it will not be longer than Nspikes/2 since we need at least two
  ## spikes to be in a burst.
  
  max.bursts <- floor(nspikes/2)
  bursts <- matrix(NA, nrow=max.bursts, ncol=3)
  colnames(bursts) = c("beg", "end", "IBI")
  burst <- 0                            #current burst number
  
  ## Phase 1 -- burst detection. Each interspike interval of the data
  ## is compared with the threshold THRE. If the interval is greater
  ## than the threshold value, it can not be part of a burst; if the
  ## interval is smaller or equal to the threhold, the interval may be
  ## part of a burst.
  
  
  
  ## LAST.END is the time of the last spike in the previous burst.
  ## This is used to calculate the IBI.
  ## For the first burst, this is no previous IBI
  last.end = NA;                        #for first burst, there is no IBI.
  
  eps<-10^(-10)
  n = 2
  in.burst = FALSE
  
  while ( n < nspikes) {
    
    next.isi = spikes[n] - spikes[n-1]
    if (in.burst) {
      if (next.isi - isi.low>eps) {
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
      if (next.isi - isi.low <=eps) {
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
  
  
  ibis = bursts[,"IBI"]
  merge.bursts = which(ibis < min.ibi)
  
  if (any(merge.bursts)) {
    ## Merge bursts efficiently.  Work backwards through the list, and
    ## then delete the merged lines afterwards.  This works when we
    ## have say 3+ consecutive bursts that merge into one.
    
    for (burst in rev(merge.bursts)) {
      bursts[burst-1, "end"] = bursts[burst, "end"]
      bursts[burst, "end"] = NA         #not needed, but helpful.
    }
    bursts = bursts[-merge.bursts,,drop=FALSE] #delete the unwanted info.
  }
  
  if (debug) {
    print("End of phase 2\n")
    print(bursts)
  }
  
  
  ## Phase 3 -- remove small bursts: less than min duration (MIN.DURN), or
  ## having too few spikes (less than MIN.SPIKES).
  ## In this phase we have the possibility of deleting all spikes.
  
  ## LEN = number of spikes in a burst.
  ## DURN = duration of burst.
  len = bursts[,"end"] - bursts[,"beg"] + 1
  durn = spikes[bursts[,"end"]] - spikes[bursts[,"beg"]]
  bursts = cbind(bursts, len, durn)
  
  rejects = which ( (durn < min.durn) | ( len < min.spikes) )
  
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