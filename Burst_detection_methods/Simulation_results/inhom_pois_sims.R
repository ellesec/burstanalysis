s<-0
u<-runif(1,0,1)
s<-s-log(u)
t<-sqrt(600*s)+10^(-16)

Lambda<-function(v){
  v-1/600*v^2
}
inhom.pois<-list()
set.seed(154)
for (i in 61:100){
s=0
v=seq(0,300,0.002)
X=0
while(X[length(X)]<=300){
  u=runif(1,0,1)
  s=s-log(u)
  t=min(v[which(Vectorize(Lambda)(v)>=s)])
  X=c(X,t)
}
X2<-X[-which(X==0)]
X2<-X2[-which(X2>300)]
isi<-diff(X2)
q<-quantile(isi, 0.1)
st<-X2[-which(isi<q)]
inhom.pois[[i]]<-st
}

i<-8
res.temp<-mclapply(inhom.pois, function(x) run.method(x, types[i], cutoffs[[i]]), mc.cores=4)

burst.stats<-sapply(inhom.pois, function(x) run.method(x, types[i], cutoffs[[i]]))

types<-c("ps", "mi", "cma", "rs", "hennig", "rgs", "logisi.pasq", "hsmm")
cutoffs<-list(-log(0.01), NULL, FALSE, -log(0.01), NULL, NULL, NULL, c(0.5, 0.9))
mi.par<-list(beg.isi=0.17,end.isi=0.3, min.ibi=0.2, min.durn=0.01, min.spikes=3)
burst.stats<-sapply(inhom.pois, function(x) run.method(x, types[1], cutoffs[[1]]))
inhom.pois.res<-data.frame(sib=unlist(burst.stats[8,]), nbursts=unlist(burst.stats[3,]), type=rep("ps",100))                
for (i in 2:8){
  if (types[i]=="rs") {
    a1<-sapply(inhom.pois, function(x) run.method(x, types[i], cutoffs[[i]]))
    sib<-do.call("rbind",a1)[,8]
    nbursts<-do.call("rbind",a1)[,3]
    inhom.pois.res<-rbind(inhom.pois.res, data.frame(sib=unlist(sib), nbursts=unlist(nbursts), type=rep(types[i], 100)))
  } else if (types[i]=="hsmm") {
    res.temp<-mclapply(inhom.pois, function(x) run.method(x, types[i], cutoffs[[i]]), mc.cores=4)
    hsmm5.sb<-sapply(res.temp, function(x) x[[1]][,8] )
    hsmm9.sb<-sapply(res.temp, function(x) x[[2]][,8] )
    hsmm5.nb<-sapply(res.temp, function(x) x[[2]][,3] )
    hsmm9.nb<-sapply(res.temp, function(x) x[[1]][,3] )
    inhom.pois.res<-rbind(inhom.pois.res, data.frame(sib=c(hsmm5.sb, hsmm9.sb), nbursts=c(hsmm5.nb, hsmm9.nb), type=c(rep("hsmm5", 100), rep("hsmm9", 100))))
    break
  } else {
    burst.stats<-sapply(inhom.pois, function(x) run.method(x, types[i], cutoffs[[i]]))
    inhom.pois.res<-rbind(inhom.pois.res, data.frame(sib=unlist(burst.stats[8,]), nbursts=unlist(burst.stats[3,]), type=rep(types[i],100)))   
  }
}

inhom.plot.subset<-subset(inhom.pois.res, type == "mi" | type == "hsmm5" | type == "logisi.pasq" | type=="ps" | type=="rs" |type=="cma")
levels(inhom.plot.subset$type)[levels(inhom.plot.subset$type)=="hsmm5"] <- "hsmm"
levels(inhom.plot.subset$type)[levels(inhom.plot.subset$type)=="logisi.pasq"] <- "logisi"
pdf("inhom_poster3.pdf", height=8, width=14)
qplot(type, sib, data = ip.sb, geom="boxplot", fill=type) + ylab("% spikes in bursts")+xlab("Burst detection method") +  scale_fill_manual(values=c(gg_color_hue(11)[c(1,4,5,7:9)]), name="", labels=c("Max Interval", "Cumulative Moving Average", "Poisson Surprise", "Rank Surprise", "logISI", "Hidden Semi-Markov Model"))+guides(fill=guide_legend(ncol=2)) +theme_bw()+theme(axis.line = element_line(colour = "black")) + theme(panel.border = element_rect(colour = "black", fill=NA, size=2))+theme(text = element_text(size=30)) + theme(legend.title= element_text(size = 30), legend.title = element_text(vjust=-1), legend.text = element_text(size = 25), legend.key.size =unit(3, "cm"), text = element_text(size=20), legend.position="right", legend.direction="vertical")
dev.off()

scale_fill_discrete(c=50, l=70, name="", labels=c("Max Interval", "Cumulative Moving Average", "Poisson Surprise", "Rank Surprise", "logISI", "Hidden Semi-Markov Model"))



pdf("inhom_pois_sib.pdf",height=8, width=15)
qplot(type, sib, data = inhomp.subset, geom="boxplot", fill=type) + ylab("% spikes in bursts")+xlab("Burst detection method") +theme_bw()+theme(axis.line = element_line(colour = "black")) + theme(panel.border = element_rect(colour = "black", fill=NA, size=2))+theme(text = element_text(size=20)) + scale_fill_discrete(c=50, l=70, name = "Burst detection method", labels=c("Poisson Surprise", "Max Interval", "Cumulative Moving Average", "Rank Surprise", "Hennig", "Robust Gaussian Surprise", "logISI", "HSMM"))+theme(legend.title= element_text(size = 20), legend.text = element_text(size = 15), legend.key.size =unit(1, "cm"), text = element_text(size=10))
dev.off()

qplot(type, sib, data = inhom.pois.res, geom="boxplot", fill=type) + ylab("% spikes in bursts")+xlab("Burst detection method") +theme_bw()+theme(axis.line = element_line(colour = "black")) + theme(panel.border = element_rect(colour = "black", fill=NA, size=2))+theme(text = element_text(size=20)) + scale_fill_discrete(c=50, l=70, name = "Burst detection method", labels=c("Poisson Surprise", "Max Interval", "Cumulative Moving Average", "Rank Surprise", "Hennig", "Robust Gaussian Surprise", "logISI", "HSMM (0.5)", "HSMM (0.9)"))+theme(legend.title= element_text(size = 20), legend.text = element_text(size = 15), legend.key.size =unit(1, "cm"), text = element_text(size=10))

inhomp.plot.subset<-subset(inhom.pois.res, type == "ps" | type == "mi" | type=="cma" | type == "logisi" | type=="hsmm5")
levels(inhomp.plot.subset$type)[levels(inhomp.plot.subset$type)=="hsmm5"] <- "hsmm"
pdf("inhom_pois_sib_pres.pdf", height=8, width=17)
qplot(type, sib, data = inhomp.plot.subset, geom="boxplot", fill=type) + ylab("% spikes in bursts")+xlab("Burst detection method") +theme_bw()+theme(axis.line = element_line(colour = "black")) + theme(panel.border = element_rect(colour = "black", fill=NA, size=2))+theme(text = element_text(size=25)) + scale_fill_discrete(c=50, l=70, name = "Burst detection method", labels=c("Poisson Surprise", "Max Interval", "Cumulative Moving Average", "logISI", "HSMM"))+theme(legend.title= element_text(size = 25), legend.text = element_text(size = 20), legend.key.size =unit(1.5, "cm"), text = element_text(size=20))
dev.off()

save(file="inhom_pois_sims.RData", inhom.pois.res, inhom.pois, inhomp.plot.subset)

pdf("inhomburst_eg.pdf", height=3, width=15)
plot(s, whichcells=1, beg=255, end=285, main="", ylab="", xlab="", xaxt="n")
axis(1, at=seq(255, 285, 5), labels=c(0,5,10,15,20,25,30))
dev.off()

plot(s, whichcells=1, beg=255, end=285, main="", ylab="", xlab="")