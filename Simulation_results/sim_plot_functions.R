#Calculates burst statistics (% spikes in bursts and number of bursts) for each burst detection method. The input is the simulated spike trains and output is data frame of burst statistics for each method
sim.result<-function(bursts){
  burst.stats<-sapply(bursts, function(x) run.method(x$spks, types[1], cutoffs[[1]]))
  res<-data.frame(sib=unlist(burst.stats[8,]), nbursts=unlist(burst.stats[3,]), type=rep("ps",100))                
  for (i in 2:length(types)){
    if (types[i]=="rs") {
      a1<-sapply(bursts, function(x) run.method(x$spks, types[i], cutoffs[[i]]))
      sib<-do.call("rbind",a1)[,8]
      nbursts<-do.call("rbind",a1)[,3]
      res<-rbind(res, data.frame(sib=unlist(sib), nbursts=unlist(nbursts), type=rep(types[i], 100)))
    } else if (types[i]=="hsmm") {
      res.temp<-mclapply(bursts, function(x) run.method(x$spks, types[i], cutoffs[[i]]), mc.cores=4)
      hsmm5.sb<-do.call("rbind",res.temp[1,])[,8]
      hsmm9.sb<-do.call("rbind",res.temp[2,])[,8]
      hsmm5.nb<-do.call("rbind",res.temp[1,])[,3]
      hsmm9.nb<-do.call("rbind",res.temp[2,])[,3]
      res<-rbind(res, data.frame(sib=c(hsmm5.sb, hsmm9.sb), nbursts=c(hsmm5.nb, hsmm9.nb), type=c(rep("hsmm5", 100), rep("hsmm9", 100))))
      break
    } else {
      burst.stats<-sapply(bursts, function(x) run.method(x$spks, types[i], cutoffs[[i]]))
      res<-rbind(res, data.frame(sib=unlist(burst.stats[8,]), nbursts=unlist(burst.stats[3,]), type=rep(types[i],100)))   
    }
  }
  res
}

#Calculates burst statistics (true positives and false positives) for each burst detection method. The input is the simulated spike trains and output is data frame of burst statistics for each method
sim.tp.result<-function(bursts){
  burst.stats<-sapply(1:100, function(x) eval.method(bursts[[x]]$spks, types[1], bursts[[x]]$burst.beg, bursts[[x]]$burst.end,  cutoffs[[1]]))
  res<-data.frame(tp=unlist(burst.stats[1,]), fp=unlist(burst.stats[2,]), type=rep("ps",100))                
  for (i in 2:length(types)){
    if (types[i]=="rs") {
      burst.stats<-sapply(1:100, function(x) eval.method(bursts[[x]]$spks, types[i], bursts[[x]]$burst.beg, bursts[[x]]$burst.end,  cutoffs[[i]]))
      tp<-do.call("rbind",burst.stats)[,1]
      fp<-do.call("rbind",burst.stats)[,2]
      res<-rbind(res, data.frame(tp=tp, fp=fp, type=rep(types[i],100)) )
    } else if (types[i]=="hsmm") {
      res.temp<-mclapply(1:100, function(x) eval.method(bursts[[x]]$spks, types[i], bursts[[x]]$burst.beg, bursts[[x]]$burst.end,  cutoffs[[i]]), mc.cores=4)
      hsmm5<-sapply(res.temp, function(x) x[[1]])
      hsmm9<-sapply(res.temp, function(x) x[[2]])
      res<-rbind(res, data.frame(tp=c(hsmm5[1,], hsmm9[1,]), fp=c(hsmm5[2,], hsmm9[2,]), type=c(rep("hsmm5", 100), rep("hsmm9", 100))))
      break
    } else {
      burst.stats<-sapply(1:100, function(x) eval.method(bursts[[x]]$spks, types[i], bursts[[x]]$burst.beg, bursts[[x]]$burst.end,  cutoffs[[i]]))
      res<-rbind(res, data.frame(tp=unlist(burst.stats[1,]), fp=unlist(burst.stats[2,]), type=rep(types[i],100)) )
    }
  }
  res
}


#ggplot options
opts<-list(theme_bw(),
           theme(axis.line = element_line(colour = "black"), axis.title.x=element_text(size=20, vjust=-0.3),axis.title.y=element_text(size=20, vjust=0.3)),
           theme(panel.border = element_rect(colour = "black", fill=NA, size=2)),
           theme(text = element_text(size=20)),
           scale_fill_discrete(c=50, l=70, name = "Burst detection method", labels=c("Poisson Surprise", "Max Interval", "Cumulative Moving Average", "Rank Surprise", "Hennig", "Robust Gaussian Surprise", "logISI", "Hidden Semi-Markov Model")),
           theme(legend.title= element_text(size = 25), legend.text = element_text(size = 20), legend.key.size =unit(2.5, "cm"), text = element_text(size=20), legend.position="none"),
           geom_hline(aes(yintercept=des.result), colour="#990000", linetype="dashed")
)

##Plot % spikes in bursts. Input is data frame of values and desired result.
plot.sib<-function(data.df, des.result) {
  data.subset<-subset(data.df, type == "ps" | type == "mi" | type=="cma" | type=="rs" | type== "hennig" |type=="rgs" | type == "logisi.pasq" | type=="hsmm5")
  levels(data.subset$type)[levels(data.subset$type)=="hsmm5"] <- "hsmm"
  levels(data.subset$type)[levels(data.subset$type)=="logisi.pasq"] <- "logisi"
  p<-qplot(type, sib, data = data.subset, geom="boxplot", fill=type) + guides(fill=guide_legend(ncol=2)) + ylab("% spikes in bursts")+xlab("Burst detection method") +opts
}

##Plot number of bursts as fraction of true number of bursts. Input is data frame of values, simulated data and desired result.
plot.nbursts<-function(data.df, bursts, des.result) {
  data.subset<-subset(data.df, type == "ps" | type == "mi" | type=="cma" | type=="rs" | type== "hennig" |type=="rgs" | type == "logisi.pasq" | type=="hsmm5")
  levels(data.subset$type)[levels(data.subset$type)=="hsmm5"] <- "hsmm"
  levels(data.subset$type)[levels(data.subset$type)=="logisi.pasq"] <- "logisi"
  num.df<-data.frame(frac.bursts=data.subset$nbursts/rep(sapply(bursts, function(x) x$num.burst),8), type=data.subset$type)
  p<-qplot(type, frac.bursts, data = num.df, geom="boxplot", fill=type) + guides(fill=guide_legend(ncol=2)) + ylab("Fraction of true number of bursts")+xlab("Burst detection method") +opts
}

##Plots fraction of true positives & false positives. Input is simulation results.
plot.tp.fp<-function(data.df) {
  data.subset<-subset(data.df, type == "ps" | type == "mi" | type=="cma" | type=="rs" | type== "hennig" |type=="rgs" | type == "logisi.pasq" | type=="hsmm5")
  levels(data.subset$type)[levels(data.subset$type)=="hsmm5"] <- "hsmm"
  levels(data.subset$type)[levels(data.subset$type)=="logisi.pasq"] <- "logisi"
  p<-list()
  p[[1]]<-qplot(type, tp, data = data.subset, geom="boxplot", fill=type) + guides(fill=guide_legend(ncol=2)) + ylab("Fraction of true positives")+xlab("Burst detection method") +opts[-7] +  geom_hline(aes(yintercept=1), colour="#990000", linetype="dashed")
  p[[2]]<-qplot(type, fp, data = data.subset, geom="boxplot", fill=type) + guides(fill=guide_legend(ncol=2)) + ylab("Fraction of false positives")+xlab("Burst detection method") +opts[-7]+  geom_hline(aes(yintercept=0), colour="#990000", linetype="dashed")
  p
}

g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

plot.sim.spikes<-function(spike.train, beg, end) {
  plot( c(beg, end), c(0,1), type='n', bty='n', yaxt="n", xaxt="n", xlab="", main="", ylab="")
  n<-length(spike.train)
  ys<-numeric(n)
  segments(spike.train, ys, spike.train, ys+0.8, lwd=0.2) 
}
