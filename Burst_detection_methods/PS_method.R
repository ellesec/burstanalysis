##PS method

PS.method<-function(spike.train, si.thresh=5) {
    si.thresh<-ifelse(is.null(si.thresh), 5, si.thresh)
    burst <- si.find.bursts.thresh(spike.train)
    if (is.null(dim(burst))){
        result<-NA
    } else {
        burst.rem<-which(burst[,"SI"]<si.thresh)
        if (length(burst.rem)) {
            burst<-burst[-burst.rem,]
        }
        if (length(dim(burst))<1) {
            burst<-data.frame(beg=burst[1], len=burst[2], SI=burst[3], durn=burst[4], mean.isis=burst[5])
        } else if (dim(burst)[1]==0){
            return(NA)
        }
        beg<-burst[,"beg"]
        len<-burst[,"len"]
        N.burst<-length(beg)
        end<-beg+len-1
        IBI<-c(NA, spike.train[beg[-1]]-spike.train[end[-N.burst]])
        result<-cbind(beg=beg, end=end, IBI=IBI, len=len, durn=burst[,"durn"], mean.isis=burst[,"mean.isis"], SI=burst[,"SI"])
        rownames(result)<-NULL
    }
    result
}




si.find.bursts.thresh<- function (spikes, debug = FALSE)
{
    nspikes = length(spikes)
    mean.isi = mean(diff(spikes))
    threshold = mean.isi/2
    n = 1
    max.bursts <- floor(nspikes/3)
    bursts <- matrix(NA, nrow = max.bursts, ncol = burst.info.len)
    burst <- 0
    while (n < nspikes - 2) {
        if (debug)
        print(n)
        if (((spikes[n + 1] - spikes[n]) < threshold) && ((spikes[n +
        2] - spikes[n + 1]) < threshold)) {
            res <- si.find.burst.thresh2(n, spikes, nspikes, mean.isi,
            burst.isi.max, debug)
            if (is.na(res[1])) {
                n <- n + 1
            }
            else {
                burst <- burst + 1
                if (burst > max.bursts) {
                    print("too many bursts")
                    browser()
                }
                bursts[burst, ] <- res
                n <- res[1] + res[2]
                names(n) <- NULL
            }
        }
        else {
            n = n + 1
        }
    }
    if (burst > 0) {
        res <- bursts[1:burst, , drop = FALSE]
        colnames(res) <- burst.info
    }
    else {
        res <- NA
    }
    res
}




si.find.burst.thresh2<-function(n, spikes, nspikes, threshold=NULL, min.si,
debug=FALSE) {
    ## Find a burst starting at spike N.
    ## Include a better phase 1.
    
    
    ## Determine ISI threshold.
    if (is.null(threshold))
    isi.thresh = 2 * mean.isi
    else
    isi.thresh = threshold
    
    if (debug)
    cat(sprintf("** find.burst %d\n", n))
    
    i=3  ## First three spikes are in burst.
    s = surprise(n, i, spikes, nspikes, mean.isi)
    
    ## Phase 1 - add spikes to the train.
    phase1 = TRUE
    ##browser()
    
    ## in Phase1, check that we still have spikes to add to the train.
    while( phase1 ) {
        
        ##printf("phase 1 s %f\n", s);
        
        i.cur = i;
        
        ## CHECK controls how many spikes we can look ahead until SI is maximised.
        ## This is normally 10, but will be less at the end of the train.
        check = min(10, nspikes-(i+n-1))
        
        looking = TRUE; okay = FALSE;
        while (looking) {
            
            if (check==0) {
                ## no more spikes left to check.
                looking=FALSE;
                break;
            }
            check=check-1; i=i+1
            s.new = surprise(n, i, spikes, nspikes, mean.isi)
            if (debug)
            printf("s.new %f s %f n %d i %d check %d\n", s.new, s, n, i, check)
            
            if (s.new > s) {
                okay=TRUE; looking=FALSE;
            } else {
                ## See if we should keep adding spikes?
                if ( (spikes[i] - spikes[i-1]) > isi.thresh ) {
                    looking = FALSE;
                }
                
            }
        }
        ## No longer checking, see if we found an improvement.
        if (okay) {
            if (s > s.new) {
                ## This should not happen.
                printf("before s %f s.new %f\n", s, s.new)
                browser()
            }
            s = s.new
        } else {
            ## Could not add more spikes onto the end of the train.
            phase1 = FALSE
            i = i.cur
        }
    }
    
    
    ## start deleting spikes from the start of the burst.
    phase2 = TRUE
    while(phase2) {
        if (i==3) {
            ## minimum length of a burst must be 3.
            phase2=FALSE
        } else {
            s.new = surprise(n+1, i-1, spikes, nspikes, mean.isi)
            if (debug)
            cat(sprintf("phase 2: n %d i %d s.new %.4f\n", n, i, s.new))
            if (s.new > s) {
                if (debug)
                print("in phase 2 acceptance\n")
                n = n+1; i = i-1
                s = s.new
            } else {
                ## removing front spike did not improve SI.
                phase2 = FALSE
            }
        }
    }
    
    
    ## End of burst detection; accumulate result.
    
    
    ## compute the ISIs, and then the mean ISI.
    
    ## Fencepost issue: I is the number of spikes in the burst, so if
    ## the first spike is N, the last spike is at N+I-1, not N+I.
    isis = diff(spikes[n+(0:(i-1))])
    mean.isis = mean(isis)
    
    durn = spikes[n+i-1] - spikes[n]
    res <- c(n=n, i=i, s=s, durn=durn, mean.isis=mean.isis)
    
    ##browser()
    res
    
}

