set.seed(129)
rate<-0.1
period<-300
mean.length<-18

long.bursts<-list()
for (i in 1:100){
  long.bursts[[i]]<-pois.burst.long(0.1, 300, 0)
}



types<-c("ps", "mi", "cma", "rs", "hennig", "rgs", "logisi.pasq", "hsmm")
cutoffs<-list(NULL, NULL, FALSE, -log(0.01), NULL, NULL, NULL, c(0.5, 0.9))
mi.par<-list(beg.isi=0.17,end.isi=0.3, min.ibi=0.2, min.durn=0.01, min.spikes=3)
burst.stats<-sapply(long.bursts, function(x) run.method(x$spks, types[1], cutoffs[[1]]))
long.burst.res<-data.frame(sib=unlist(burst.stats[8,]), nbursts=unlist(burst.stats[3,]), type=rep("ps",100))                
for (i in 2:8){
  if (types[i]=="rs") {
    a1<-sapply(long.bursts, function(x) run.method(x$spks, types[i], cutoffs[[i]]))
    sib<-do.call("rbind",a1)[,8]
    nbursts<-do.call("rbind",a1)[,3]
    long.burst.res<-rbind(long.burst.res, data.frame(sib=unlist(sib), nbursts=unlist(nbursts), type=rep(types[i], 100)))
  } else if (types[i]=="hsmm") {
    res.temp<-mclapply(long.bursts, function(x) run.method(x$spks, types[i], cutoffs[[i]]), mc.cores=4)
    hsmm5.sb<-sapply(res.temp, function(x) x[[1]][,8] )
    hsmm9.sb<-sapply(res.temp, function(x) x[[2]][,8] )
    hsmm5.nb<-sapply(res.temp, function(x) x[[2]][,3] )
    hsmm9.nb<-sapply(res.temp, function(x) x[[1]][,3] )
   long.burst.res<-rbind(long.burst.res, data.frame(sib=c(hsmm5.sb, hsmm9.sb), nbursts=c(hsmm5.nb, hsmm9.nb), type=c(rep("hsmm5", 100), rep("hsmm9", 100))))
    break
  } else {
    burst.stats<-sapply(long.bursts, function(x) run.method(x$spks, types[i], cutoffs[[i]]))
    long.burst.res<-rbind(long.burst.res, data.frame(sib=unlist(burst.stats[8,]), nbursts=unlist(burst.stats[3,]), type=rep(types[i],100)))   
  }
}

lburst.num.df<-data.frame(frac.bursts=long.burst.res$nbursts/rep(sapply(long.bursts, function(x) x$num.burst),10), type=long.burst.res$type)
#PLOTS
qplot(type, sib, data = long.burst.res, geom="boxplot", fill=type) + ylab("% spikes in bursts")+xlab("Burst detection method") + theme(legend.position="none")+theme_bw()+theme(axis.line = element_line(colour = "black"))+ scale_fill_discrete(c=50, l=70, labels=c("Poisson Surprise", "Max Interval", "Cumulate Moving Average", "Rank Surprise", "Hennig", "Robust Gaussian Surprise", "logISI", "HSMM (0.5)", "HSMM (0.9)"))
qplot(type, frac.bursts, data = lburst.num.df, geom="boxplot", fill=type) + ylab("Fraction of true bursts")+xlab("Burst detection method") + theme(legend.position="none")+theme_bw()+theme(axis.line = element_line(colour = "black"))+ scale_fill_discrete(c=50, l=70, labels=c("Poisson Surprise", "Max Interval", "Cumulate Moving Average", "Rank Surprise", "Hennig", "Robust Gaussian Surprise", "logISI", "HSMM (0.5)", "HSMM (0.9)"))

long.burst.subset<-subset(long.burst.res, type == "ps" | type == "mi" | type=="cma" | type=="rs" | type== "hennig" |type=="rgs" | type == "logisi.pasq" | type=="hsmm5")
levels(long.burst.subset$type)[levels(long.burst.subset$type)=="hsmm5"] <- "hsmm"
levels(long.burst.subset$type)[levels(long.burst.subset$type)=="logisi.pasq"] <- "logisi"

###TO PLOT
pdf("long_burst_sib.pdf",height=8, width=15)
qplot(type, sib, data = long.burst.res, geom="boxplot", fill=type) + ylab("% spikes in bursts")+xlab("Burst detection method") +theme_bw()+theme(axis.line = element_line(colour = "black")) + theme(panel.border = element_rect(colour = "black", fill=NA, size=2))+theme(text = element_text(size=20)) + scale_fill_discrete(c=50, l=70, name = "Burst detection method", labels=c("Poisson Surprise", "Max Interval", "Cumulative Moving Average", "Rank Surprise", "Hennig", "Robust Gaussian Surprise", "logISI", "HSMM (0.5)", "HSMM (0.9)"))+theme(legend.title= element_text(size = 20), legend.text = element_text(size = 15), legend.key.size =unit(1, "cm"), text = element_text(size=10))
dev.off()

lburst.num.df<-data.frame(frac.bursts=long.burst.res$nbursts/rep(sapply(long.bursts, function(x) x$num.burst),10), type=long.burst.res$type)
pdf("long_burst_fracb.pdf",height=8, width=15)
qplot(type, frac.bursts, data = lburst.num.df, geom="boxplot", fill=type) + ylab("Fraction of true bursts")+xlab("Burst detection method") +theme_bw()+theme(axis.line = element_line(colour = "black")) + theme(panel.border = element_rect(colour = "black", fill=NA, size=2))+theme(text = element_text(size=20)) + scale_fill_discrete(c=50, l=70, name = "Burst detection method", labels=c("Poisson Surprise", "Max Interval", "Cumulative Moving Average", "Rank Surprise", "Hennig", "Robust Gaussian Surprise", "logISI", "HSMM (0.5)", "HSMM (0.9)"))+theme(legend.title= element_text(size = 20), legend.text = element_text(size = 15), legend.key.size =unit(1, "cm"), text = element_text(size=10))
dev.off()

pdf("longburst_eg.pdf", height=3, width=15)
s$spikes[[1]]<-long.bursts[[1]]$spks
plot(s, whichcells=1, beg=0, end=50, main="", ylab="", xlab="")
dev.off()

###PRESENTATION PLOTS
longburst.plot.subset<-subset(long.burst.res, type == "ps" | type == "mi" | type=="cma" | type == "logisi" | type=="hsmm5")
levels(longburst.plot.subset$type)[levels(longburst.plot.subset$type)=="hsmm5"] <- "hsmm"
pdf("longburst_sib_pres.pdf", height=8, width=17)
qplot(type, sib, data = longburst.plot.subset, geom="boxplot", fill=type) + ylab("% spikes in bursts")+xlab("Burst detection method") +theme_bw()+theme(axis.line = element_line(colour = "black")) + theme(panel.border = element_rect(colour = "black", fill=NA, size=2))+theme(text = element_text(size=25)) + scale_fill_discrete(c=50, l=70, name = "Burst detection method", labels=c("Poisson Surprise", "Max Interval", "Cumulative Moving Average", "logISI", "HSMM"))+theme(legend.title= element_text(size = 25), legend.text = element_text(size = 20), legend.key.size =unit(1.5, "cm"), text = element_text(size=20))
dev.off()

lburst.plot<-subset(lburst.num.df, type == "ps" | type == "mi" | type=="cma" | type == "logisi" | type=="hsmm5")
levels(lburst.plot$type)[levels(lburst.plot$type)=="hsmm5"] <- "hsmm"
pdf("longburst_nb_pres.pdf", height=8, width=17)
qplot(type, frac.bursts, data = lburst.plot, geom="boxplot", fill=type)  +geom_hline(yintercept=1, linetype="dashed", color = "black", size=1)+ylab("Fraction of true bursts detected")+xlab("Burst detection method") +theme_bw()+theme(axis.line = element_line(colour = "black")) + theme(panel.border = element_rect(colour = "black", fill=NA, size=2))+theme(text = element_text(size=25)) + scale_fill_discrete(c=50, l=70, name = "Burst detection method", labels=c("Poisson Surprise", "Max Interval", "Cumulative Moving Average", "logISI", "HSMM"))+theme(legend.title= element_text(size = 25), legend.text = element_text(size = 20), legend.key.size =unit(1.5, "cm"), text = element_text(size=20))
dev.off()




##PLOTS for presentation
longburst.plot.subset<-subset(long.burst.res, type == "ps" | type == "mi" | type=="cma" | type == "logisi" | type=="hsmm5")
levels(longburst.plot.subset$type)[levels(longburst.plot.subset$type)=="hsmm5"] <- "hsmm"


pdf("long_burst_sib.pdf", height=8, width=15)
qplot(type, sib, data = long.burst.res, geom="boxplot", fill=type) + ylab("% spikes in bursts")+xlab("Burst detection method") +theme_bw()+theme(axis.line = element_line(colour = "black")) + theme(panel.border = element_rect(colour = "black", fill=NA, size=2))+theme(text = element_text(size=20)) + scale_fill_discrete(c=50, l=70, name = "Burst detection method", labels=c("Poisson Surprise", "Max Interval", "Cumulate Moving Average", "Rank Surprise", "Hennig", "Robust Gaussian Surprise", "logISI", "HSMM (0.5)", "HSMM (0.9)")+theme(legend.title= element_text(size = 20), legend.text = element_text(size = 20), legend.key.size =unit(1.5, "cm"), text = element_text(size=15))
dev.off()


pdf("longburst_eg.pdf", height=3, width=15)
plot(s, whichcells=1, beg=0, end=60, main="", ylab="", xlab="")
dev.off()

save(file="long_burst_sims.RData", long.bursts, long.burst.res, lburst.num.df)

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

pdf("rast.pdf", height=20, width=30)
par(mfrow=c(3,6))
for (i in 1:18){
  plot(s.list[[18+i]])
}



