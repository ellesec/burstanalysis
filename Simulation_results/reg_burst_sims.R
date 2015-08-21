set.seed(128)
rate<-0.2
period<-300
mean.length<-5

reg.bursts<-list()
for (i in 1:100){
  reg.bursts[[i]]<-pois.burst.noiseless(0.2, 300, 0)
}

h<-hist(log(diff(reg.bursts[[1]]$spks)),40)
logisi.method(reg.bursts[[2]]$spks)

logisi.compute.train(reg.bursts[[1]]$spks)$Locmin

types<-c("ps", "mi", "cma", "rs", "hennig", "rgs", "logisi.pasq", "hsmm")
cutoffs<-list(NULL, NULL, FALSE, -log(0.01), NULL, NULL, NULL, c(0.5, 0.9))
mi.par<-list(beg.isi=0.17,end.isi=0.3, min.ibi=0.2, min.durn=0.01, min.spikes=3)
burst.stats<-sapply(reg.bursts, function(x) run.method(x$spks, types[1], cutoffs[[1]]))
reg.burst.res<-data.frame(sib=unlist(burst.stats[8,]), nbursts=unlist(burst.stats[3,]), type=rep("ps",100))                
for (i in 2:8){
  if (types[i]=="rs") {
    a1<-sapply(reg.bursts, function(x) run.method(x$spks, types[i], cutoffs[[i]]))
    sib<-do.call("rbind",a1)[,8]
    nbursts<-do.call("rbind",a1)[,3]
    reg.burst.res<-rbind(reg.burst.res, data.frame(sib=unlist(sib), nbursts=unlist(nbursts), type=rep(types[i], 100)))
  } else if (types[i]=="hsmm") {
    res.temp<-sapply(reg.bursts, function(x) run.method(x$spks, types[i], cutoffs[[i]]))
    hsmm5.sb<-do.call("rbind",res.temp[1,])[,8]
    hsmm9.sb<-do.call("rbind",res.temp[2,])[,8]
    hsmm5.nb<-do.call("rbind",res.temp[1,])[,3]
    hsmm9.nb<-do.call("rbind",res.temp[2,])[,3]
    reg.burst.res<-rbind(reg.burst.res, data.frame(sib=c(hsmm5.sb, hsmm9.sb), nbursts=c(hsmm5.nb, hsmm9.nb), type=c(rep("hsmm5", 100), rep("hsmm9", 100))))
    break
  } else {
    burst.stats<-sapply(reg.bursts, function(x) run.method(x$spks, types[i], cutoffs[[i]]))
    reg.burst.res<-rbind(reg.burst.res, data.frame(sib=unlist(burst.stats[8,]), nbursts=unlist(burst.stats[3,]), type=rep(types[i],100)))   
  }
}

save(file="reg_burst_sims.RData", reg.burst.res, reg.bursts)
#PLOTS

reg.plot.subset<-subset(reg.burst.res, type == "mi" | type == "hsmm5" | type == "logisi.pasq" | type=="ps" | type=="rs" |type=="cma")
levels(reg.plot.subset$type)[levels(reg.plot.subset$type)=="hsmm5"] <- "hsmm"
levels(reg.plot.subset$type)[levels(reg.plot.subset$type)=="logisi.pasq"] <- "logisi"
pdf("reg_burst_poster2.pdf", height=8, width=14)
qplot(type, sib, data = rb.ps, geom="boxplot", fill=type) + ylab("% spikes in bursts")+xlab("Burst detection method") +theme_bw()+theme(axis.line = element_line(colour = "black")) + theme(panel.border = element_rect(colour = "black", fill=NA, size=2))+theme(text = element_text(size=30)) + scale_fill_discrete(c=50, l=70, name = "Burst detection method", labels=c("Poisson Surprise", "Max Interval", "Cumulative Moving Average", "logISI", "HSMM"))+theme(legend.title= element_text(size = 25), legend.text = element_text(size = 20), legend.key.size =unit(1.5, "cm"), text = element_text(size=30),  legend.position="none")
dev.off()



pdf("reg_burst_sib.pdf", height=8, width=15)
qplot(type, sib, data = reg.burst.res, geom="boxplot", fill=type) + ylab("% spikes in bursts")+xlab("Burst detection method") + theme(legend.position="none")+theme_bw()+theme(axis.line = element_line(colour = "black"))+ scale_fill_discrete(c=50, l=70, name = "Burst detection method", labels=c("Poisson Surprise", "Max Interval", "Cumulate Moving Average", "Rank Surprise", "Hennig", "Robust Gaussian Surprise", "logISI", "HSMM (0.5)", "HSMM (0.9)")) +  theme(panel.border = element_rect(colour = "black", fill=NA, size=2))+theme(text = element_text(size=20)) + theme(legend.title= element_text(size = 16), legend.text = element_text(size = 15), legend.key.size =unit(1, "cm"), text = element_text(size=15))
dev.off()

burst.num.df<-data.frame(frac.bursts=reg.burst.res$nbursts/rep(sapply(reg.bursts, function(x) x$num.burst),10), type=reg.burst.res$type)
qplot(type, frac.bursts, data = burst.num.df, geom="boxplot", fill=type) + ylab("Fraction of true bursts")+xlab("Burst detection method") + theme(legend.position="none")+theme_bw()+theme(axis.line = element_line(colour = "black"))+ scale_fill_discrete(c=50, l=70, labels=c("Poisson Surprise", "Max Interval", "Cumulate Moving Average", "Rank Surprise", "Hennig", "Robust Gaussian Surprise", "logISI", "HSMM (0.5)", "HSMM (0.9)"))


###PRESENTATION PLOTS
regburst.plot.subset<-subset(reg.burst.res, type == "ps" | type == "mi" | type=="cma" | type == "logisi" | type=="hsmm5")
levels(regburst.plot.subset$type)[levels(regburst.plot.subset$type)=="hsmm5"] <- "hsmm"
pdf("regburst_pres.pdf", height=8, width=17)
qplot(type, sib, data = regburst.plot.subset, geom="boxplot", fill=type) + ylab("% spikes in bursts")+xlab("Burst detection method") +theme_bw()+theme(axis.line = element_line(colour = "black")) + theme(panel.border = element_rect(colour = "black", fill=NA, size=2))+theme(text = element_text(size=25)) + scale_fill_discrete(c=50, l=70, name = "Burst detection method", labels=c("Poisson Surprise", "Max Interval", "Cumulative Moving Average", "logISI", "HSMM"))+theme(legend.title= element_text(size = 25), legend.text = element_text(size = 20), legend.key.size =unit(1.5, "cm"), text = element_text(size=20))
dev.off()



pdf("regburst_eg.pdf", height=3, width=15)
s$spikes[[1]]<-reg.bursts[[8]]$spks
plot(s, whichcells=1, beg=0, end=30, main="", ylab="", xlab="")
dev.off()

pois.burst.noiseless<-function(rate, period, start) {
  burst.times<-sort(runif(rpois(1,period*rate), min=start, max=start+period))
  no.bursts<-length(burst.times)
   burst.length<-rpois(no.bursts, mean.length)
  burst.pos<-sapply(burst.length, function(x) runif(x, min=-0.15, max=0.15))
  burst.pos.sorted<-lapply( burst.pos, sort)
  spike.times<-list()
  for ( i in 1:length(burst.times)){
    spike.times[[i]]<-burst.pos.sorted[[i]]+burst.times[i]
  }
  spike.times2<-spike.times[which(burst.length>2)]
  burst.beg<-sapply(spike.times2, min)
  burst.end<-sapply(spike.times2, max)
  overlap<-burst.beg[-1]-burst.end[-length(burst.end)]
  ol<-which(overlap<0.5)
  if (length(ol)>0) {
    spike.times2<-spike.times2[-ol]
  }
  burst.beg<-sapply(spike.times2, min)
  burst.end<-sapply(spike.times2, max)
  sts<-sort(unlist(spike.times2))
  sts<-sts[which(sts>0)]
  
  insim<-list()
  insim$spks<-sts
  insim$burst.beg<-burst.beg
  insim$burst.end<-burst.end
  insim$num.bursts<-length(burst.beg)
  insim
}

##PLOT
pdf("reg_burst_sib.pdf", height=8, width=15)
qplot(type, sib, data = reg.burst.res, geom="boxplot", fill=type) + ylab("% spikes in bursts")+xlab("Burst detection method") +theme_bw()+theme(axis.line = element_line(colour = "black")) + theme(panel.border = element_rect(colour = "black", fill=NA, size=2))+theme(text = element_text(size=20)) + scale_fill_discrete(c=50, l=70, name = "Burst detection method", labels=c("Poisson Surprise", "Max Interval", "Cumulative Moving Average", "Rank Surprise", "Hennig", "Robust Gaussian Surprise", "logISI", "HSMM (0.5)", "HSMM (0.9)"))+theme(legend.title= element_text(size = 20), legend.text = element_text(size = 15), legend.key.size =unit(1, "cm"), text = element_text(size=10))
dev.off()

pdf("regburst_eg.pdf", height=3, width=15)
s$spikes[[1]]<-reg.bursts[[1]]$spks
plot(s, whichcells=1, beg=0, end=30, main="", ylab="", xlab="")
dev.off()




