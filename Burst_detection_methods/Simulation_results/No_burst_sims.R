set.seed(123)
period<-300
rate<-0.5
start<-0
end<-0
pois.no.bursts2<-list()
for ( i in 1:50) {
  spike.times<-sort(runif(rpois(1,period*rate), min=start, max=start+period))
  isi<-diff(spike.times)
  rem.spikes<-which(isi<quantile(isi, 0.1))
pois.no.bursts[[i]]<-spike.times[-rem.spikes]
}


set.seed(126)
period<-300
gamma.no.bursts<-list()
for (i in 1:50) {
  gst<-rgamma(200, 1,  rate=0.5)
  q<-quantile(gst, 0.1)
  gst<-gst[-which(gst<q)]
  gsum<-cumsum(gst)
  gsum<-gsum[which(gsum<period)]
gamma.no.bursts[[i]]<-gsum
}

all.no.bursts<-append(pois.no.bursts, gamma.no.bursts)

types<-c("ps", "mi", "cma", "rs", "hennig", "rgs", "logisi.pasq", "hsmm")
cutoffs<-list(-log(0.01), NULL, FALSE, -log(0.01), NULL, NULL, NULL, c(0.5, 0.9))
mi.par<-list(beg.isi=0.17,end.isi=0.3, min.ibi=0.2, min.durn=0.01, min.spikes=3)

no.burst.res<-data.frame(unlist(sapply(all.no.bursts, function(x) run.method(x, types[1], cutoffs[[1]]))[8,]))
names(no.burst.res)[1]<-types[1]                 
for (i in 2:8){
  if (types[i]=="rs") {
    a1<-sapply(all.no.bursts, function(x) run.method(x, types[i], cutoffs[[i]]))
  a<-do.call("rbind",a1)[,8]
  no.burst.res<-cbind(no.burst.res, unlist(a))
  names(no.burst.res)[i]<-types[i]
} else if (types[i]=="hsmm") {
  res.temp<-sapply(all.no.bursts, function(x) run.method(x, types[i], cutoffs[[i]]))
  hsmm5<-do.call("rbind",res.temp[1,])[,8]
  hsmm9<-do.call("rbind",res.temp[2,])[,8]
  no.burst.res<-cbind(no.burst.res, hsmm5=hsmm5, hsmm9=hsmm9)
  break
} else {
a<-unlist(sapply(all.no.bursts, function(x) run.method(x, types[i], cutoffs[[i]]))[8,])
no.burst.res<-cbind(no.burst.res, unlist(a))
names(no.burst.res)[i]<-types[i]
}
}

boxplot(no.burst.res, col=Set2)

no.burst.res2<-data.frame(val=unlist(sapply(all.no.bursts, function(x) run.method(x, types[1], cutoffs[[1]]))[8,]), type=rep("ps",100))                
for (i in 2:8){
  if (types[i]=="rs") {
    a1<-sapply(all.no.bursts, function(x) run.method(x, types[i], cutoffs[[i]]))
    a<-do.call("rbind",a1)[,8]
    no.burst.res2<-rbind(no.burst.res2, data.frame(val=unlist(a), type=rep(types[i], 100)))
  } else if (types[i]=="hsmm") {
    res.temp<-sapply(all.no.bursts, function(x) run.method(x, types[i], cutoffs[[i]]))
    hsmm5<-do.call("rbind",res.temp[1,])[,8]
    hsmm9<-do.call("rbind",res.temp[2,])[,8]
    no.burst.res2<-rbind(no.burst.res2, data.frame(val=c(hsmm5, hsmm9), type=c(rep("hsmm5", 100), rep("hsmm9", 100))))
    break
  } else {
    a<-unlist(sapply(all.no.bursts, function(x) run.method(x, types[i], cutoffs[[i]]))[8,])
    no.burst.res2<-rbind(no.burst.res2, data.frame(val=unlist(a), type=rep(types[i], 100)))
  }
}

#plots
no.burst.subset<-subset(no.burst.res2, type == "ps" | type == "mi" | type=="cma" | type=="rs" | type== "hennig" |type=="rgs" | type == "logisi.pasq" | type=="hsmm5")
levels(no.burst.subset$type)[levels(no.burst.subset$type)=="hsmm5"] <- "hsmm"
levels(no.burst.subset$type)[levels(no.burst.subset$type)=="logisi.pasq"] <- "logisi"
qplot(type, val, data = no.burst.res2, geom="boxplot", fill=type) + ylab("% spikes in bursts")+xlab("Burst detection method") + theme(legend.position="none")+theme_bw()+theme(axis.line = element_line(colour = "black"))+ scale_fill_discrete(c=50, l=70, labels=c("Poisson Surprise", "Max Interval", "Cumulative Moving Average", "Rank Surprise", "Hennig", "Robust Gaussian Surprise", "logISI", "HSMM (0.5)", "HSMM (0.9)"))
pdf("no_burst_sub.pdf", height=8, width=15)
qplot(type, val, data = no.burst.subset, geom="boxplot", fill=type) + ylab("% spikes in bursts")+xlab("Burst detection method") +theme_bw()+theme(axis.line = element_line(colour = "black")) + theme(panel.border = element_rect(colour = "black", fill=NA, size=2))+theme(text = element_text(size=25)) + scale_fill_discrete(c=50, l=70, name = "Burst detection method", labels=c("Poisson Surprise", "Max Interval", "Cumulative Moving Average", "Rank Surprise", "Hennig", "Robust Gaussian Surprise", "logISI", "HSMM"))+theme(legend.title= element_text(size = 20), legend.text = element_text(size = 15), legend.key.size =unit(1, "cm"), text = element_text(size=20))
dev.off()


##Presentation plots
noburst.plot.subset<-subset(no.burst.res2, type == "mi" | type == "hsmm5" | type == "logisi.pasq" | type=="ps" | type=="rs" |type=="cma")
levels(noburst.plot.subset$type)[levels(noburst.plot.subset$type)=="hsmm5"] <- "hsmm"
levels(noburst.plot.subset$type)[levels(noburst.plot.subset$type)=="logisi.pasq"] <- "logisi"
pdf("nb_poster3.pdf", height=8, width=14)
qplot(type, val, data = nb.ps3, geom="boxplot", fill=type) + ylab("Percentage of spikes in bursts")+xlab("Burst detection method") +theme_bw()+theme(axis.line = element_line(colour = "black")) + theme(panel.border = element_rect(colour = "black", fill=NA, size=2))+theme(text = element_text(size=30)) + scale_fill_manual(values=c(gg_color_hue(11)[c(1,4,5,7:9)]), name="", labels=c("Max Interval", "Cumulative Moving Average", "Poisson Surprise", "Rank Surprise", "logISI", "Hidden Semi-Markov Model"))+theme(legend.title= element_text(size = 25), legend.text = element_text(size = 20), legend.key.size =unit(1.5, "cm"), text = element_text(size=30),  legend.position="none")
dev.off()

position_jitter( height = NULL)

pdf("noburst_eg.pdf", height=3, width=15)
s$spikes[[1]]<-all.no.bursts[[2]]
plot(s, whichcells=1, beg=0, end=30, main="", ylab="", xlab="")
dev.off()



save(file="no_burst_sims.RData", no.burst.res2, all.no.bursts, no.burst.ray, no.burst.res3)

no.burst.ray<-list()
set.seed(324)
for (i in 1:100) {
isi.sim<-rgenray(200, 3, 1)
st<-cumsum(isi.sim)
no.burst.ray[[i]]<-st[which(st<300)]
}

no.burst.res3<-data.frame(val=unlist(sapply(no.burst.ray, function(x) run.method(x, types[1], cutoffs[[1]]))[8,]), type=rep("ps",100))                
for (i in 2:8){
  if (types[i]=="rs") {
    a1<-sapply(no.burst.ray, function(x) run.method(x, types[i], cutoffs[[i]]))
    a<-do.call("rbind",a1)[,8]
    no.burst.res3<-rbind(no.burst.res3, data.frame(val=unlist(a), type=rep(types[i], 100)))
  } else if (types[i]=="hsmm") {
    res.temp<-sapply(no.burst.ray, function(x) run.method(x, types[i], cutoffs[[i]]))
    hsmm5<-do.call("rbind",res.temp[1,])[,8]
    hsmm9<-do.call("rbind",res.temp[2,])[,8]
    no.burst.res3<-rbind(no.burst.res3, data.frame(val=c(hsmm5, hsmm9), type=c(rep("hsmm5", 100), rep("hsmm9", 100))))
    break
  } else {
    a<-unlist(sapply(no.burst.ray, function(x) run.method(x, types[i], cutoffs[[i]]))[8,])
    no.burst.res3<-rbind(no.burst.res3, data.frame(val=unlist(a), type=rep(types[i], 100)))
  }
}

pdf("no_burst_ray.pdf", height=8, width=15)
qplot(type, val, data = no.burst.res3, geom="boxplot", fill=type) + ylab("% spikes in bursts")+xlab("Burst detection method") +theme_bw()+theme(axis.line = element_line(colour = "black")) + theme(panel.border = element_rect(colour = "black", fill=NA, size=2))+theme(text = element_text(size=20)) + scale_fill_discrete(c=50, l=70, name = "Burst detection method", labels=c("Poisson Surprise", "Max Interval", "Cumulative Moving Average", "Rank Surprise", "Hennig", "Robust Gaussian Surprise", "logISI", "HSMM (0.5)", "HSMM (0.9)"))+theme(legend.title= element_text(size = 20), legend.text = element_text(size = 15), legend.key.size =unit(1, "cm"), text = element_text(size=10))
dev.off()