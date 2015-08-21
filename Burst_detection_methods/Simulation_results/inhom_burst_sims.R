set.seed(135)
rate.seq<-seq(0.3, 0.8, 0.1)
mean.lengths<-20*rate.seq

inhom.bursts<-list()
for (i in 1:100){
total<-list()
N<-0
for (j in 1:6) {
  a<-pois.burst.inhom(rate.seq[j], 50, 50*(j-1), mean.lengths[j], rate.seq[j])
  total$spks<-c(total$spks, a$spks)
  total$burst.beg<-c(total$burst.beg, (a$burst.beg+N))
  total$burst.end<-c(total$burst.end, (a$burst.end+N) )
  N<-N+length(a$spks)
}
inhom.bursts[[i]]<-total
}

types<-c("ps", "mi", "cma", "rs", "hennig", "rgs", "logisi.pasq", "hsmm")
cutoffs<-list(NULL, NULL, FALSE, -log(0.01), NULL, NULL, NULL, c(0.5, 0.9))
mi.par<-list(beg.isi=0.17,end.isi=0.3, min.ibi=0.2, min.durn=0.01, min.spikes=3)
burst.stats<-sapply(1:100, function(x) eval.method(inhom.bursts[[x]]$spks, types[1], inhom.bursts[[x]]$burst.beg, inhom.bursts[[x]]$burst.end,  cutoffs[[1]]))
inhom.burst.res<-data.frame(tp=unlist(burst.stats[1,]), fp=unlist(burst.stats[2,]), type=rep("ps",100))                
for (i in 2:7){
  if (types[i]=="rs") {
    burst.stats<-sapply(1:100, function(x) eval.method(inhom.bursts[[x]]$spks, types[i], inhom.bursts[[x]]$burst.beg, inhom.bursts[[x]]$burst.end,  cutoffs[[i]]))
    tp<-do.call("rbind",burst.stats)[,1]
    fp<-do.call("rbind",burst.stats)[,2]
    inhom.burst.res<-rbind(inhom.burst.res, data.frame(tp=tp, fp=fp, type=rep(types[i],100)) )
  } else if (types[i]=="hsmm") {
    res.temp<-mclapply(1:100, function(x) eval.method(inhom.bursts[[x]]$spks, types[i], inhom.bursts[[x]]$burst.beg, inhom.bursts[[x]]$burst.end,  cutoffs[[i]]), mc.cores=4)
    hsmm5<-sapply(res.temp, function(x) x[[1]])
    hsmm9<-sapply(res.temp, function(x) x[[2]])
    inhom.burst.res<-rbind(inhom.burst.res, data.frame(tp=c(hsmm5[1,], hsmm9[1,]), fp=c(hsmm5[2,], hsmm9[2,]), type=c(rep("hsmm5", 100), rep("hsmm9", 100))))
    break
  } else {
    burst.stats<-sapply(1:100, function(x) eval.method(inhom.bursts[[x]]$spks, types[i], inhom.bursts[[x]]$burst.beg, inhom.bursts[[x]]$burst.end,  cutoffs[[i]]))
    inhom.burst.res<-rbind(inhom.burst.res, data.frame(tp=unlist(burst.stats[1,]), fp=unlist(burst.stats[2,]), type=rep(types[i],100)) )
  }
}

qplot(type, tp, data = inhom.burst.res, geom="boxplot", fill=type)

pois.burst.inhom<-function(rate, period, start, mean.length, noise.rate ) {
  burst.times<-sort(runif(rpois(1,period*rate), min=start, max=start+period))
  no.bursts<-length(burst.times)
  #burst.length<-NULL
  #for (i in 1:no.bursts) {
  #  mn<- round(burst.times[i], -2)/100+4
  #  burst.length[i]<-rpois(1, mn)
  #}
  burst.length<-rpois(no.bursts, mean.length)
  burst.pos<-sapply(burst.length, function(x) runif(x, min=-0.5, max=0.5))
  burst.pos.sorted<-lapply( burst.pos, sort)
  spike.times<-list()
  for ( i in 1:length(burst.times)){
    spike.times[[i]]<-burst.pos.sorted[[i]]+burst.times[i]
  }
  spike.times2<-spike.times[which(burst.length>2)]
  while (spike.times2[[1]][1]<start){
    spike.times2<-spike.times2[-1]
  } 
  N.temp<-length(spike.times2)
  while (spike.times2[[N.temp]][length(spike.times2[[N.temp]])]>(period+start)) {
    spike.times2<-spike.times2[-N.temp]
    N.temp<-length(spike.times2)
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

  
  
  noise.times<-sort(runif(rpois(1,period*noise.rate), min=start, max=start+period))
  isi<-diff(noise.times)
  rem.spikes<-which(isi<quantile(isi, 0.1))
  noise.filtered<-noise.times[-rem.spikes]
  
  
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

pdf("inhom_burst_tp.pdf", height=8, width=15)
qplot(type, tp, data = inhom.burst.res, geom="boxplot", fill=type) + ylab("Fraction of true positives")+xlab("Burst detection method") +theme_bw()+theme(axis.line = element_line(colour = "black")) + theme(panel.border = element_rect(colour = "black", fill=NA, size=2))+theme(text = element_text(size=20)) + scale_fill_discrete(c=50, l=70, name = "Burst detection method", labels=c("Poisson Surprise", "Max Interval", "Cumulative Moving Average", "Rank Surprise", "Hennig", "Robust Gaussian Surprise", "logISI", "HSMM (0.5)", "HSMM (0.9)"))+theme(legend.title= element_text(size = 20), legend.text = element_text(size = 15), legend.key.size =unit(1, "cm"), text = element_text(size=10))
dev.off()

pdf("inhom_burst_fp.pdf", height=8, width=15)
qplot(type, fp, data = inhom.burst.res, geom="boxplot", fill=type) + ylab("Fraction of false positives")+xlab("Burst detection method") +theme_bw()+theme(axis.line = element_line(colour = "black")) + theme(panel.border = element_rect(colour = "black", fill=NA, size=2))+theme(text = element_text(size=20)) + scale_fill_discrete(c=50, l=70, name = "Burst detection method", labels=c("Poisson Surprise", "Max Interval", "Cumulative Moving Average", "Rank Surprise", "Hennig", "Robust Gaussian Surprise", "logISI", "HSMM (0.5)", "HSMM (0.9)"))+theme(legend.title= element_text(size = 20), legend.text = element_text(size = 15), legend.key.size =unit(1, "cm"), text = element_text(size=10))
dev.off()

pdf("inhomburst_eg.pdf", height=3, width=15)
s$spikes[[1]]<-inhom.bursts[[1]]$spks
plot(s, whichcells=1, beg=0, end=200, main="", ylab="", xlab="")
dev.off()


save(file="inhom_burst_sims.RData", inhom.bursts, inhom.burst.res)

