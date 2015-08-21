HSMM.method <- function(spike.train, cutoff, scan = c(5,25), dur.shape = c(1, 1)){
  
  # Load and execute burstHSMM package
  require(burstHSMM)
  y <- diff(spike.train)
  if (length(y)<2) {
      Results<-rep(list(NA), length(cutoff))
      return(Results)
  }
  a <- fit.hsmm(y, scan, dur.shape, limit = NULL,
                states = NULL, pars = NULL, hpar = NULL, tune = NULL,
                nburn = 0, nsamp = 1e3, nskip = 1, nshot = 5)
  
  # Extract final probabilities
  A <- summary(a,toplot=F)
  p <- A$prob.burst
  
  # Determine results for all cutoff values
  Results <- list()
  for (x in 1:length(cutoff)){
    below.cutoff <- ifelse(p <= cutoff[x],1,0)
    spike.ID <- diff(below.cutoff)
    n <- length(below.cutoff)
    if (below.cutoff[1] == 1){
      start.pos <- which(spike.ID==-1) +1
    }else{
      start.pos <- c(1, which(spike.ID==-1) +1)
    }
    if (below.cutoff[n] == 1){
      end.pos <- which(spike.ID==1) +1
    }else{
      end.pos <- c(which(spike.ID==1) +1, n)
    }
    N.burst<-length(start.pos)
    if (N.burst<1) {
      result<-NA
    } else{
    start.times<-spike.train[start.pos]
    end.times<-spike.train[end.pos]
    beg<-start.pos
    end<-end.pos
    IBI<-c(NA, start.times[-1]-end.times[-N.burst])
    len<-end.pos-start.pos+1
    durn<-end.times-start.times
    mean.isis<-durn/(len-1)
    result<-cbind(beg, end, IBI, len, durn, mean.isis, SI=rep(cutoff[x], N.burst))
    }
    Results[[x]] <- result
  }
  return(Results)
  
}