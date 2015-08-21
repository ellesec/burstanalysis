#Bursts with rate 0.5, length 8, on [-0.4, 0.4]. Noise gamma(0.5, 1)
set.seed(130)
rate<-0.5
period<-300
mean.length<-8

noisy.bursts<-list()
for (i in 1:100){
  noisy.bursts[[i]]<-pois.burst.noisy(0.5,300,0)
}

types<-c("ps", "mi", "cma", "rs", "hennig", "rgs", "logisi.pasq", "hsmm")
cutoffs<-list(NULL, NULL, FALSE, -log(0.01), NULL, NULL, NULL, c(0.5, 0.9))
mi.par<-list(beg.isi=0.17,end.isi=0.3, min.ibi=0.2, min.durn=0.01, min.spikes=3)
burst.stats<-sapply(1:100, function(x) eval.method(noisy.bursts[[x]]$spks, types[1], noisy.bursts[[x]]$burst.beg, noisy.bursts[[x]]$burst.end,  cutoffs[[1]]))
noisy.burst.res<-data.frame(tp=unlist(burst.stats[1,]), fp=unlist(burst.stats[2,]), type=rep("ps",100))                
for (i in 2:7){
  if (types[i]=="rs") {
    burst.stats<-sapply(1:100, function(x) eval.method(noisy.bursts[[x]]$spks, types[i], noisy.bursts[[x]]$burst.beg, noisy.bursts[[x]]$burst.end,  cutoffs[[i]]))
    tp<-do.call("rbind",burst.stats)[,1]
    fp<-do.call("rbind",burst.stats)[,2]
    noisy.burst.res<-rbind(noisy.burst.res, data.frame(tp=tp, fp=fp, type=rep(types[i],100)) )
  } else if (types[i]=="hsmm") {
    res.temp<-mclapply(1:100, function(x) eval.method(noisy.bursts[[x]]$spks, types[i], noisy.bursts[[x]]$burst.beg, noisy.bursts[[x]]$burst.end,  cutoffs[[i]]), mc.cores=4)
    hsmm5<-sapply(res.temp, function(x) x[[1]])
    hsmm9<-sapply(res.temp, function(x) x[[2]])
    noisy.burst.res<-rbind(noisy.burst.res, data.frame(tp=c(hsmm5[1,], hsmm9[1,]), fp=c(hsmm5[2,], hsmm9[2,]), type=c(rep("hsmm5", 100), rep("hsmm9", 100))))
    break
  } else {
    burst.stats<-sapply(1:100, function(x) eval.method(noisy.bursts[[x]]$spks, types[i], noisy.bursts[[x]]$burst.beg, noisy.bursts[[x]]$burst.end,  cutoffs[[i]]))
    noisy.burst.res<-rbind(noisy.burst.res, data.frame(tp=unlist(burst.stats[1,]), fp=unlist(burst.stats[2,]), type=rep(types[i],100)) )
  }
}

noisy.burst.res<-rbind(noisy.burst.res, nbr)
save(file="noisy_burst_sims.RData", noisy.burst.res, noisy.bursts, res.temp)

##To plot

pdf("noisy_burst_tp.pdf", height=8, width=15)
qplot(type, tp, data = noisy.burst.res, geom="boxplot", fill=type) + ylab("Fraction of true positives")+xlab("Burst detection method") +theme_bw()+theme(axis.line = element_line(colour = "black")) + theme(panel.border = element_rect(colour = "black", fill=NA, size=2))+theme(text = element_text(size=20)) + scale_fill_discrete(c=50, l=70, name = "Burst detection method", labels=c("Poisson Surprise", "Max Interval", "Cumulative Moving Average", "Rank Surprise", "Hennig", "Robust Gaussian Surprise", "logISI", "HSMM (0.5)", "HSMM (0.9)"))+theme(legend.title= element_text(size = 20), legend.text = element_text(size = 15), legend.key.size =unit(1, "cm"), text = element_text(size=10))
dev.off()

pdf("noisy_burst_fp.pdf", height=8, width=15)
qplot(type, fp, data = noisy.burst.res, geom="boxplot", fill=type) + ylab("Fraction of false positives")+xlab("Burst detection method") +theme_bw()+theme(axis.line = element_line(colour = "black")) + theme(panel.border = element_rect(colour = "black", fill=NA, size=2))+theme(text = element_text(size=20)) + scale_fill_discrete(c=50, l=70, name = "Burst detection method", labels=c("Poisson Surprise", "Max Interval", "Cumulative Moving Average", "Rank Surprise", "Hennig", "Robust Gaussian Surprise", "logISI", "HSMM (0.5)", "HSMM (0.9)"))+theme(legend.title= element_text(size = 20), legend.text = element_text(size = 15), legend.key.size =unit(1, "cm"), text = element_text(size=10))
dev.off()

pdf("noisyburst_eg.pdf", height=3, width=15)
s$spikes[[1]]<-noisy.bursts[[5]]$spks
plot(s, whichcells=1, beg=0, end=30, main="", ylab="", xlab="")
dev.off()

##PRESENTATION PLOTS
noisyburst.plot.subset<-subset(noisy.burst.res, type == "ps" | type == "mi" | type=="cma" | type == "logisi" | type=="hsmm5")
levels(noisyburst.plot.subset$type)[levels(noisyburst.plot.subset$type)=="hsmm5"] <- "hsmm"

pdf("noisyburst_tp_pres.pdf", height=8, width=17)
qplot(type, tp, data = noisyburst.plot.subset, geom="boxplot", fill=type) + ylab("Fraction of true positives")+xlab("Burst detection method") +theme_bw()+theme(axis.line = element_line(colour = "black")) + theme(panel.border = element_rect(colour = "black", fill=NA, size=2))+theme(text = element_text(size=25)) + scale_fill_discrete(c=50, l=70, name = "Burst detection method", labels=c("Poisson Surprise", "Max Interval", "Cumulative Moving Average", "logISI", "HSMM"))+theme(legend.title= element_text(size = 25), legend.text = element_text(size = 20), legend.key.size =unit(1.5, "cm"), text = element_text(size=20))
dev.off()  

pdf("noisyburst_fp_pres.pdf", height=8, width=17)
qplot(type, fp, data = noisyburst.plot.subset, geom="boxplot", fill=type) + ylab("Fraction of false positives")+xlab("Burst detection method") +theme_bw()+theme(axis.line = element_line(colour = "black")) + theme(panel.border = element_rect(colour = "black", fill=NA, size=2))+theme(text = element_text(size=25)) + scale_fill_discrete(c=50, l=70, name = "Burst detection method", labels=c("Poisson Surprise", "Max Interval", "Cumulative Moving Average", "logISI", "HSMM"))+theme(legend.title= element_text(size = 25), legend.text = element_text(size = 20), legend.key.size =unit(1.5, "cm"), text = element_text(size=20))
dev.off() 




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
