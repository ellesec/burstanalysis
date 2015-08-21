
gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=45, c=80)[1:n]
}

plot(PS.ROC9[,3], PS.ROC9[,2], col=colors[1], "l", xlab="1-Specificity", ylab="Sensitivity", xlim=c(0,0.1), ylim=c(0,1), lwd=2)
lines(MI.ROC9[,3], MI.ROC9[,2], col=colors[2], lwd=2, lty=2)
points(CMA.ROC9[1,3], CMA.ROC9[1,2], pch=16, col=colors[3])
points(log.ROC9[,3], log.ROC9[,2], pch=16, col=colors[4])
legend("bottomright",c("PS", "MI", "CMA", "logisi"),  bty="n", col=c(colors[1], colors[2], colors[3], colors[4]), lty = c(1,1,1,1,1,1,1,1))
lines(HSMM.ROC9[,3], HSMM.ROC9[,2], col=colors[5], lwd=2)

pdf("ROC_curves.pdf", width=14, height=7) 
par(mfrow=c(1,2))
plot(PS.ROC11[,3], PS.ROC11[,2], col=colors[1], "l", xlab="1-Specificity", ylab="Sensitivity", xlim=c(0,0.9), ylim=c(0,1), lwd=3)
lines(MI.ROC11[,3], MI.ROC11[,2], col=colors[2], lwd=3)
points(CMA.ROC11[1,3], CMA.ROC11[1,2], pch=16, col=colors[3])
points(log.ROC11[,3], log.ROC11[,2], pch=16, col=colors[4])
lines(HSMM.ROC11[,3], HSMM.ROC11[,2],col=colors[5], lwd=3)
legend("bottomright",c("PS", "MI", "CMA", "logisi", "HSMM"),  bty="n", col=c(colors[1], colors[2], colors[3], colors[4], colors[5]), lty = c(1,1,1,1,1,1,1,1))

ROC9<-cbind(PS.ROC9[,2:3], type="ps")
ROC9<-rbind(ROC9, data.frame(tp=MI.ROC9[,2],fp=MI.ROC9[,3], type="mi"))
ROC9<-rbind(ROC9, data.frame(tp=CMA.ROC9[1,2],fp=CMA.ROC9[1,3], type="cma"))
ROC9<-rbind(ROC9, data.frame(tp=log.ROC9[2],fp=log.ROC9[3], type="logisi"))
ROC9<-rbind(ROC9, cbind(HSMM.ROC9[,2:3], type="hsmm"))  
dummy<-p11 +theme(legend.position=c(0.2 ,0.2), legend.direction="vertical")
dummy.legend <- g_legend(dummy)


pdf("ROCs.pdf", height=6, width=20)
grid.arrange(p11, p13, ncol=3)
pushViewport(viewport(0.9, 0.6, 0.2, 0.4))
grid.draw(dummy.legend); popViewport()
pushViewport(viewport())
grid.text(LETTERS[1:2], gp=gpar(fontsize=20),
          x=unit(c(0.02,  0.35), "npc"),
          y=unit(c(0.98, 0.98), "npc"))
dev.off()

p_samp<-ggplot(data=ROC11[,], aes(x=fp, y=tp, group=type, color=type))+geom_line()  + scale_fill_manual(values=colour.hue.gg(8,70,50)[c(1, 2,4,6, 8)], name = "Burst detection method", labels=c("Poisson Surprise", "Max Interval", "Cumulative Moving Average", "logISI", "HSMM"))+theme(legend.title= element_text(size = 25), legend.text = element_text(size = 20), legend.key.size =unit(2.5, "cm"), text = element_text(size=20), legend.position="none")


p11<-ggplot(data=ROC11[,], aes(x=fp, y=tp, group=type, color=type))+geom_line(size=1.3) +theme_bw()+geom_point(data = ROC11[121,], colour = colour.hue.gg(5, 70, 50)[3], size=3)+ geom_point(data = ROC11[122,], colour = colour.hue.gg(5,70,50)[4], size=3) + ylab("Sensitivity")+xlab("1-Specificity") + theme(axis.line = element_line(colour = "black"), axis.title.x=element_text(size=20, vjust=-0.3),axis.title.y=element_text(size=20, vjust=0.3)) + theme(panel.border = element_rect(colour = "black", fill=NA, size=2))+theme(text = element_text(size=20)) + scale_colour_discrete( c=50, l=70, name  ="Burst detection method",                                                                    labels=c("Poisson Surprise", "Max Interval", "Cumulative Moving Average", "logISI", "HSMM"))+theme(legend.title= element_text(size = 25), legend.text = element_text(size = 20), legend.key.size =unit(2.5, "cm"), text = element_text(size=20), legend.position="none")
p13<-ggplot(data=ROC13[-c(138, 139, 140, 142, 143),], aes(x=fp, y=tp, group=type, color=type))+geom_line(size=1.2) +theme_bw()+geom_point(data = ROC13[124,], colour = colour.hue.gg(5, 70, 50)[3], size=3)+ geom_point(data = ROC13[125,], colour = colour.hue.gg(5,70,50)[4], size=3) + ylab("Sensitivity")+xlab("1-Specificity") + theme(axis.line = element_line(colour = "black"), axis.title.x=element_text(size=20, vjust=-0.3),axis.title.y=element_text(size=20, vjust=0.3)) + theme(panel.border = element_rect(colour = "black", fill=NA, size=2))+theme(text = element_text(size=20)) + scale_colour_discrete( c=50, l=70, name  ="Payer",                                                                    labels=c("Poisson Surprise", "Max Interval", "Cumulative Moving Average", "logISI", "HSMM"))+theme(legend.title= element_text(size = 25), legend.text = element_text(size = 20), legend.key.size =unit(2.5, "cm"), text = element_text(size=20), legend.position="none")




p11.paper<-ggplot(data=ROC11[,], aes(x=fp, y=tp, group=type, color=type))+geom_line(size=1.3) +theme_bw()+geom_point(data = ROC11[121,], colour = colour.hue.gg(8, 70, 50)[3], size=3)+ geom_point(data = ROC11[122,], colour = colour.hue.gg(8,70,50)[7], size=3) + ylab("Sensitivity")+xlab("1-Specificity") + theme(axis.line = element_line(colour = "black"), axis.title.x=element_text(size=20, vjust=-0.3),axis.title.y=element_text(size=20, vjust=0.3)) + theme(panel.border = element_rect(colour = "black", fill=NA, size=2))+theme(text = element_text(size=20)) + scale_color_manual(values=colour.hue.gg(8,70,50)[c(1, 2,3,7,8)])+theme(legend.title= element_text(size = 25), legend.text = element_text(size = 20), legend.key.size =unit(2.5, "cm"), text = element_text(size=20), legend.position="none")

p13.paper<-ggplot(data=ROC13[-c(138, 139, 140, 142, 143),], aes(x=fp, y=tp, group=type, color=type))+geom_line(size=1.3) +theme_bw()+geom_point(data = ROC13[124,], colour = colour.hue.gg(8, 70, 50)[3], size=3)+ geom_point(data = ROC13[125,], colour = colour.hue.gg(8,70,50)[7], size=3) + ylab("Sensitivity")+xlab("1-Specificity") + theme(axis.line = element_line(colour = "black"), axis.title.x=element_text(size=20, vjust=-0.3),axis.title.y=element_text(size=20, vjust=0.3)) + theme(panel.border = element_rect(colour = "black", fill=NA, size=2))+theme(text = element_text(size=20)) + scale_color_manual(values=colour.hue.gg(8,70,50)[c(1, 2,3,7,8)])+theme(legend.title= element_text(size = 25), legend.text = element_text(size = 20), legend.key.size =unit(2.5, "cm"), text = element_text(size=20), legend.position="none")

p15.paper<-ggplot(data=ROC15[-c(7,9,11,12,13, 14, 15,26),], aes(x=fp, y=tp, group=type, color=type))+geom_line(size=1.3) +theme_bw()+geom_point(data = ROC15[100,], colour = colour.hue.gg(8, 70, 50)[3], size=3)+ geom_point(data = ROC15[101,], colour = colour.hue.gg(8,70,50)[7], size=3) + ylab("Sensitivity")+xlab("1-Specificity") + theme(axis.line = element_line(colour = "black"), axis.title.x=element_text(size=20, vjust=-0.3),axis.title.y=element_text(size=20, vjust=0.3)) + theme(panel.border = element_rect(colour = "black", fill=NA, size=2))+theme(text = element_text(size=20)) + scale_color_manual(values=colour.hue.gg(8,70,50)[c(1, 2,3,7,8)])+theme(legend.title= element_text(size = 25), legend.text = element_text(size = 20), legend.key.size =unit(2.5, "cm"), text = element_text(size=20), legend.position="none")

pdf("paper_ROC.pdf", height=6, width=20)
grid.arrange(p11.paper, p13.paper, p15.paper, ncol=3)
grid.text(LETTERS[1:3], gp=gpar(fontsize=20),
          x=unit(c(0.02,  0.35, 0.68), "npc"),
          y=unit(c(0.98, 0.98, 0.98), "npc"))
dev.off()


scale_fill_manual(values=colour.hue.gg(8,70,50)[c(1, 2,3,7)], name = "Burst detection method", labels=c("Max Interval", "Poisson Surprise", "Cumulative Moving Average", "logISI"))

ggplot(data=ROC15[-c(7,9,11,12,13, 14, 15,26),], aes(x=fp, y=tp, group=type, color=type))+geom_line(size=1.2) +theme_bw()+geom_point(data = ROC15[100,], colour = colour.hue.gg(5, 70, 50)[3], size=3)+ geom_point(data = ROC15[101,], colour = colour.hue.gg(5,70,50)[4], size=3) + ylab("Sensitivity")+xlab("1-Specificity") + theme(axis.line = element_line(colour = "black"), axis.title.x=element_text(size=20, vjust=-0.3),axis.title.y=element_text(size=20, vjust=0.3)) + theme(panel.border = element_rect(colour = "black", fill=NA, size=2))+theme(text = element_text(size=20)) + scale_colour_discrete( c=50, l=70, name  ="Payer",                                                                    labels=c("Poisson Surprise", "Max Interval", "Cumulative Moving Average", "logISI", "HSMM"))+theme(legend.title= element_text(size = 25), legend.text = element_text(size = 20), legend.key.size =unit(2.5, "cm"), text = element_text(size=20), legend.position="none")

plot(PS.ROC13[,3], PS.ROC13[,2], col=colors[1], "l", xlab="1-Specificity", ylab="Sensitivity", xlim=c(0,0.9), ylim=c(0,1), lwd=3)
lines(MI.ROC13[,2], MI.ROC13[,1], col=colors[2], lwd=3)
points(CMA.ROC13[1,3], CMA.ROC13[1,2], pch=16, col=colors[3])
points(LOG.ROC13[3], LOG.ROC13[2], pch=16, col=colors[4])
lines(HSMM.ROC13[-c(13:16, 18:19),3], HSMM.ROC13[-c(13:16, 18:19),2],col=colors[5], lwd=3)
legend("bottomright",c("PS", "MI", "CMA", "logisi", "HSMM"),  bty="n", col=c(colors[1], colors[2], colors[3], colors[4], colors[5]), lty = c(1,1,1,1,1,1,1,1))

text(grconvertX(c(0.02, 0.52), from='ndc'),
     grconvertY(c(0.9, 0.9), from='ndc'),
     c('A',  'B'), xpd=NA, cex=1.5, font=2)
dev.off()

for (ms in seq(0.2, 3, 0.05)){
  mi.par$beg.isi<-ms
  mi.par$end.isi<-ms+0.13
  mi.res<-rbind( mi.res, c(apply(sapply(s11.indx, function(x) eval.method(s11$spikes[[x]], "mi", s11$allb[[x]][,1], s11$allb[[x]][,2], i)),1, mean)))
}

for ( i in 1:length(p.vals)){
result<-eval.stats(spike.train, beg.true, end.true, bursts)
}  

plot(PS.ROC11[,3], PS.ROC11[,2], col=colors[1], "l", xlab="1-Specificity", ylab="Sensitivity", xlim=c(0,0.15), ylim=c(0,1), lwd=2)
d<-cbind(pvals, do.call("rbind",a))

s13.indx<-c(8, 15, 20, 27, 29)
a<-sapply(s13.indx, function(x) eval.method(s13$spikes[[x]], "logisi.pasq", s13$allb[[x]][,1], s13$allb[[x]][,2], TRUE))
CMA.ROC13<-rbind(CMA.ROC13, c(FALSE, apply(a, 1, mean)))
cbind(ps.res, sapply(s13.indx, function(x) eval.method(s13$spikes[[x]], "ps", s13$allb[[x]][,1], s13$allb[[x]][,2], i))
for (i in seq(5, 200, 5)) {
  ps.res<-cbind(ps.res, apply(sapply(s13.indx, function(x) eval.method(s13$spikes[[x]], "ps", s13$allb[[x]][,1], s13$allb[[x]][,2], i)),1, mean))
}
lines(HSMM.ROC11[,3], HSMM.ROC11[,2], col=colors[5], lwd=2)



0.737069 0.000000
0.8814815 0.0000000
0.9653580 0.1428571
0.651255 0.000000
0.8869048 0.0000000


si.bursts<-spikes.to.bursts.ed(s11, "si")
s9.indx<-c(8,11,16,22,24)
s11.indx<-c(1,2,3,5,6)
s15.indx<-c(17,12,26,24)

thresh<-10
PS.ROC15<-NULL
for (j in seq(10, 250, 10)){
  res<-data.frame()
  for (i in s15.indx) {
    st<-s15$spikes[[i]]
    ps.res<-PS.method.thresh(st, j)
    #ps.res<-cbind(beg=si.bursts[[i]][,1], end=si.bursts[[i]][,1]+si.bursts[[i]][,2]-1, SI=si.bursts[[i]][,3] )
    #indx<-which(ps.res[,"SI"]>j)
    #ps.res<-cbind(beg=ps.res[indx,1], end=ps.res[indx,2])
    if (!is.na(ps.res[1])) {
   
    ev<-eval.stats(st, s15$allb[[i]][,1],s15$allb[[i]][,2], ps.res )
    } else {
      ev<-c(0,0)
    }
    res<-rbind(res, ev)
  }
  
  PS.ROC15<-rbind(PS.ROC15, data.frame(val=j, tp= mean(res[,1], na.rm=TRUE), fp=mean(res[,2], na.rm=TRUE)))
}

RGS.ROC9<-NULL
for (j in seq(0.05,1, 0.055)){
  res<-data.frame()
  for (i in s9.indx[3:5]) {
    st<-s9$spikes[[i]]
    ps.res<-RGS.method(st, j)
    #ps.res<-cbind(beg=si.bursts[[i]][,1], end=si.bursts[[i]][,1]+si.bursts[[i]][,2]-1, SI=si.bursts[[i]][,3] )
    #indx<-which(ps.res[,"SI"]>j)
    #ps.res<-cbind(beg=ps.res[indx,1], end=ps.res[indx,2])
    if (!is.na(ps.res[1])) {
      ev<-eval.stats(st, s9$allb[[i]][,1],s9$allb[[i]][,2], ps.res )
    } else {
      ev<-c(0,0)
    }
    res<-rbind(res, ev)
  }
  
  RGS.ROC9<-rbind(RGS.ROC9, data.frame(val=j, tp= mean(res[,1], na.rm=TRUE), fp=mean(res[,2], na.rm=TRUE)))
}







for (cutoff in seq(0.8, 0.90, 0.1)){
  
  res<-data.frame()
  for (i in 1:39) {
    st<-s15$spikes[[i]]
    isi.list[[1]]<-data.frame(isi = c(st[1], diff(st)))
    res1<-f.BPsummary(data=isi.list,Pthresh=cutoff)$burst[[1]][,7:8]
    start.true<- st[s15$allb[[i]][,"beg"]]
    end.true<- st[s15$allb[[i]][,"end"]]
    
    #start.method<- RS.res[[i]][[1]][,1]
    #end.method<- RS.res[[i]][[1]][,2]
    start.method<-res1[,1]
    end.method<-res1[,2]
    res<-rbind(res, evaluate(st, start.true, end.true, start.method, end.method))
  }
  
  RGS.ROC15<-rbind(RGS.ROC15, data.frame(val=cutoff, tp= mean(res[,1], na.rm=TRUE), fp=mean(res[,2], na.rm=TRUE)))
  
}

for (cutoff in seq(0.01, 0.09, 0.01)){
  RS.res<-lapply(s15$spikes, function(x) ranksurprisemethod(x, -log(cutoff)))
  res<-data.frame()
  for (i in 1:39) {
    st<-s15$spikes[[i]]
    isi.list[[1]]<-data.frame(isi = c(st[1], diff(st)))
    start.true<- st[s15$allb[[i]][,"beg"]]
    end.true<- st[s15$allb[[i]][,"end"]]
    
    start.method<- RS.res[[i]][[1]][,1]
    end.method<- RS.res[[i]][[1]][,2]
    #start.method<-res1[,1]
    #end.method<-res1[,2]
    res<-rbind(res, eval.stats(st, s9$allb[[i]][,1],s9$allb[[i]][,2], ps.res ))
  }
  
  RS.ROC15<-rbind(RS.ROC15, data.frame(val=cutoff, tp= mean(res[,1], na.rm=TRUE), fp=mean(res[,2], na.rm=TRUE)))
  
}

0.00001 0.1176367 0.09965244
32 0.00005 0.1465458 0.25086212
31 0.00010 0.1485169 0.25086212
30 0.00050 0.1888897 0.25086212
29 0.00100 0.1925401 0.25086212
28 0.00500

for (j in 1:100) {
  st<-9$spikes[[j]]
  isi.list[[1]]<-data.frame(isi = c(st[1], diff(st)))
  cutoff<-NULL
  res<-f.BPsummary(data=isi.list,Pthresh=cutoff)$burst[[1]][,7:8]
  RGS.tb[[j]]<-res
}

hb<-lapply(s15$spikes, hennig.burst)

#for (j in seq(5, 200, 5)){
res<-data.frame()
for (i in 1:39) {
  st<-s15$spikes[[i]]
  start.true<- st[s15$allb[[i]][,"beg"]]
  end.true<- st[s15$allb[[i]][,"end"]]
  #si.adj<-subset(si.bursts[[i]], si.bursts[[i]][,"SI"]>j)
  #si.adj<-si.bursts[[i]]
  start.method<- hb[[i]][,1]
  end.method<- hb[[i]][,2]
  res<-rbind(res, evaluate(st, start.true, end.true, start.method, end.method))
}

hen.ROC15<-rbind(hen.ROC15, data.frame(val=0.15, tp= mean(res[,1]), fp=mean(res[,2], na.rm=TRUE)))


}

for (pv in seq(0.3, 0.9, 0.1)) {
  
  pvals<-c(0.01, 0.025, 0.05, seq(0.1, 0.9, 0.1), seq(0.91, 1, 0.01))
  HSb15<-mclapply(s15$spikes[s15.indx], function(x) hiddensemimarkovmodelmethod(x, pvals), mc.cores=4)
  for (j in 1 :21){
    #for (j in seq(5, 200, 5)){
    res<-data.frame()
    for (ii in 1:4) {
      i<-s15.indx[[ii]]
      st<-s15$spikes[[i]]
      #si.adj<-subset(si.bursts[[i]], si.bursts[[i]][,"SI"]>j)
      #si.adj<-si.bursts[[i]]
      hs.res<-HSb15[[ii]][[j]]
      if (hs.res[1]>0){
        ps.res<-cbind(beg = sapply(hs.res[,1], function(x) which(st==x)), end = sapply(hs.res[,2], function(x) which(st==x)))
      res<-rbind(res, eval.stats(st, s15$allb[[i]][,1],s15$allb[[i]][,2], ps.res ))
      }
      else {
        res<-rbind(res, c(0,0))
      }
    }
    
    HSMM.ROC15<-rbind(HSMM.ROC15, data.frame(val=pvals[j], tp= mean(res[,1], na.rm=TRUE), fp=mean(res[,2], na.rm=TRUE)))
    
  }
  
  res<-data.frame()
  cmab<-list()
  isi.low <- logisi.compute(s15)$Locmin
  logisi.par$isi.low <- 0.02
  ncells <- s$NCells
  for (i in 1:39){
    st<-s15$spikes[[i]]
    start.true<- st[s15$allb[[i]][,"beg"]]
    end.true<- st[s15$allb[[i]][,"end"]]
    logb[[i]]<-logisi.find.burst(s15$spikes[[i]])
    #cmab[[i]]<-CMA.bursts(s15$spikes[[i]])
    if (is.null(dim(logb[[i]]))){
      start.method<- NULL
      end.method<- NULL
    } else {
      start.method<- st[logb[[i]][,1]]
      end.method<- st[logb[[i]][,2]]
    }
    res<-rbind(res, evaluate(st, start.true, end.true, start.method, end.method))
  }
  
  log.ROC15<-data.frame(val=NA, tp= mean(res[-13,1], na.tm=TRUE), fp=mean(res[-13,2], na.rm=TRUE))
  
 MI.ROC9<-NULL
  for (ms in seq(0.2, 3, 0.05)){
    mi.par$beg.isi<-ms
    mi.par$end.isi<-ms+0.13
    res<-data.frame()
    for (i in s9.indx) {
      st<-s9$spikes[[i]]
      #si.adj<-subset(si.bursts[[i]], si.bursts[[i]][,"SI"]>j)
      mi.res<-MI.method(st)
      res<-rbind(res, eval.stats(st,s9$allb[[i]][,"beg"], s9$allb[[i]][,"end"], mi.res ))
    }
    
    
    MI.ROC9<-rbind(MI.ROC9, data.frame(val=ms, tp= mean(res[,1]), fp=mean(res[,2], na.rm=TRUE)))
    
  }
 
 log.ROC9<-NULL
 res<-data.frame()
 for (i in s9.indx) {
   st<-s9$spikes[[i]]
   #si.adj<-subset(si.bursts[[i]], si.bursts[[i]][,"SI"]>j)
   mi.res<-logisi.method(st)
  #mi.res<-log.res[[i]]
   res<-rbind(res, eval.stats(st,s9$allb[[i]][,"beg"], s9$allb[[i]][,"end"], mi.res ))
 }
log.ROC9<-rbind(log.ROC9, data.frame(val=TRUE, tp= mean(res[,1]), fp=mean(res[,2], na.rm=TRUE)))
  

CMA.ROC9<-NULL
res<-data.frame()
for (i in s9.indx) {
  st<-s9$spikes[[i]]
  #si.adj<-subset(si.bursts[[i]], si.bursts[[i]][,"SI"]>j)
  mi.res<-CMA.method(st, TRUE)
  #mi.res<-log.res[[i]]
  res<-rbind(res, eval.stats(st,s9$allb[[i]][,"beg"], s9$allb[[i]][,"end"], mi.res ))
}
CMA.ROC9<-rbind(CMA.ROC9, data.frame(val=TRUE, tp= mean(res[,1]), fp=mean(res[,2], na.rm=TRUE)))

 
  pdf("ROC_Demas9.pdf", width=10, height=10 )
  plot(PS.ROC[,3], PS.ROC[,2], "l", xlab="1-Specificity", ylab="Sensitivity", col=p[1], xlim=c(0,0.5), ylim=c(0,1), main = "Demas2003_P9")
  lines(MI.ROC[,3], MI.ROC[,2], col=p[3])
  points(CMA.ROC[,3], CMA.ROC[,2], col=p[6], pch=16)
  lines(hen.ROC[,3], hen.ROC[,2], col=p[7])
  points(log.ROC[,3], log.ROC[,2], col=p[9], pch=16)
  lines(RGS.ROC[,3], RGS.ROC[,2], col=p[11])
  lines(RS.ROC[,3], RS.ROC[,2], col=p[13])
  lines(HSMM.ROC[7:18,3], HSMM.ROC[7:18,2], col=p[15])
  legend("bottomright",c("PS", "MI", "CMA", "hen", "logisi", "RGS", "RS", "HSMM"),  bty="n", col=c(p[1], p[3], p[6], p[7], p[9], p[11], p[13], p[15]), lty = c(1,1,1,1,1,1,1,1))
  dev.off()
  
  pdf("ROC_Demas11.pdf", width=10, height=10 )
  plot(PS.ROC11[,3], PS.ROC11[,2], "l", xlab="1-Specificity", ylab="Sensitivity", col=p[1], xlim=c(0,0.6), ylim=c(0,1), main = "Demas2003_P11")
  lines(MI.ROC11[,3], MI.ROC11[,2], col=p[3])
  points(CMA.ROC11[,3], CMA.ROC11[,2], col=p[6], pch=16)
  lines(hen.ROC11[,3], hen.ROC11[,2], col=p[7])
  points(log.ROC11[,3], log.ROC11[,2], col=p[9], pch=16)
  lines(RGS.ROC11[,3], RGS.ROC11[,2], col=p[11])
  lines(RS.ROC11[,3], RS.ROC11[,2], col=p[13])
  lines(HSMM.ROC11[7:18,3], HSMM.ROC11[7:18,2], col=p[15])
  legend("bottomright",c("PS", "MI", "CMA", "hen", "logisi", "RGS", "RS", "HSMM"),  bty="n", col=c(p[1], p[3], p[6], p[7], p[9], p[11], p[13], p[15]), lty = c(1,1,1,1,1,1,1,1))
  dev.off()
  
  pdf("ROC_Demas13.pdf", width=10, height=10 )
  plot(PS.ROC13[,3], PS.ROC13[,2], "l", xlab="1-Specificity", ylab="Sensitivity", col=p[1], xlim=c(0,0.4), ylim=c(0,1), main = "Demas2003_P13")
  lines(MI.ROC13[,3], MI.ROC13[,2], col=p[3])
  points(CMA.ROC13[,3], CMA.ROC13[,2], col=p[6], pch=16)
  lines(hen.ROC13[,3], hen.ROC13[,2], col=p[7])
  points(log.ROC13[,3], log.ROC13[,2], col=p[9], pch=16)
  lines(RGS.ROC13[,3], RGS.ROC13[,2], col=p[11])
  lines(RS.ROC13[,3], RS.ROC13[,2], col=p[13])
  lines(HSMM.ROC13[7:18,3], HSMM.ROC13[7:18,2], col=p[15])
  legend("bottomright",c("PS", "MI", "CMA", "hen", "logisi", "RGS", "RS", "HSMM"),  bty="n", col=c(p[1], p[3], p[6], p[7], p[9], p[11], p[13], p[15]), lty = c(1,1,1,1,1,1,1,1))
  dev.off()
  
  pdf("ROC_Demas15.pdf", width=10, height=10 )
  plot(PS.ROC15[,3], PS.ROC15[,2], "l", xlab="1-Specificity", ylab="Sensitivity", col=p[1], xlim=c(0,0.4), ylim=c(0,1), main = "Demas2003_P15")
  lines(MI.ROC15[,3], MI.ROC15[,2], col=p[3])
  points(CMA.ROC15[,3], CMA.ROC15[,2], col=p[6], pch=16)
  lines(hen.ROC15[,3], hen.ROC15[,2], col=p[7])
  points(log.ROC15[,3], log.ROC15[,2], col=p[9], pch=16)
  lines(RGS.ROC15[,3], RGS.ROC15[,2], col=p[11])
  lines(RS.ROC15[,3], RS.ROC15[,2], col=p[13])
  lines(HSMM.ROC15[,3], HSMM.ROC15[,2], col=p[15])
  legend("bottomright",c("PS", "MI", "CMA", "hen", "logisi", "RGS", "RS", "HSMM"),  bty="n", col=c(p[1], p[3], p[6], p[7], p[9], p[11], p[13], p[15]), lty = c(1,1,1,1,1,1,1,1))
  dev.off()
  
  pdf("ROC_Demas.pdf", width=12, height=15 )
  par(mfrow=c(2,2), las=1)
  plot(PS.ROC[,3], PS.ROC[,2], "l", xlab="1-Specificity", ylab="Sensitivity", col=p[1], xlim=c(0,0.6), ylim=c(0,1), main = "")
  lines(MI.ROC[,3], MI.ROC[,2], col=p[3])
  points(CMA.ROC[,3], CMA.ROC[,2], col=p[6], pch=16)
  lines(hen.ROC[,3], hen.ROC[,2], col=p[7])
  points(log.ROC[,3], log.ROC[,2], col=p[9], pch=16)
  lines(RGS.ROC[,3], RGS.ROC[,2], col=p[11])
  lines(RS.ROC[,3], RS.ROC[,2], col=p[13])
  lines(HSMM.ROC[6:18,3], HSMM.ROC[6:18,2], col=p[15])
  legend("bottomright",c("PS", "MI", "CMA", "hen", "logisi", "RGS", "RS", "HSMM"),  bty="n", col=c(p[1], p[3], p[6], p[7], p[9], p[11], p[13], p[15]), lty = c(1,1,1,1,1,1,1,1))
  plot(PS.ROC11[,3], PS.ROC11[,2], "l", xlab="1-Specificity", ylab="Sensitivity", col=p[1], xlim=c(0,0.6), ylim=c(0,1), main = "")
  lines(MI.ROC11[,3], MI.ROC11[,2], col=p[3])
  points(CMA.ROC11[,3], CMA.ROC11[,2], col=p[6], pch=16)
  lines(hen.ROC11[,3], hen.ROC11[,2], col=p[7])
  points(log.ROC11[,3], log.ROC11[,2], col=p[9], pch=16)
  lines(RGS.ROC11[,3], RGS.ROC11[,2], col=p[11])
  lines(RS.ROC11[,3], RS.ROC11[,2], col=p[13])
  lines(HSMM.ROC11[,3], HSMM.ROC11[,2], col=p[15])
  legend("bottomright",c("PS", "MI", "CMA", "hen", "logisi", "RGS", "RS", "HSMM"),  bty="n", col=c(p[1], p[3], p[6], p[7], p[9], p[11], p[13], p[15]), lty = c(1,1,1,1,1,1,1,1))
  plot(PS.ROC13[,3], PS.ROC13[,2], "l", xlab="1-Specificity", ylab="Sensitivity", col=p[1], xlim=c(0,0.6), ylim=c(0,1), main = "")
  lines(MI.ROC13[,3], MI.ROC13[,2], col=p[3])
  points(CMA.ROC13[,3], CMA.ROC13[,2], col=p[6], pch=16)
  lines(hen.ROC13[,3], hen.ROC13[,2], col=p[7])
  points(log.ROC13[,3], log.ROC13[,2], col=p[9], pch=16)
  lines(RGS.ROC13[,3], RGS.ROC13[,2], col=p[11])
  lines(RS.ROC13[,3], RS.ROC13[,2], col=p[13])
  lines(HSMM.ROC13[3:19,3], HSMM.ROC13[3:19,2], col=p[15])
  legend("bottomright",c("PS", "MI", "CMA", "hen", "logisi", "RGS", "RS", "HSMM"),  bty="n", col=c(p[1], p[3], p[6], p[7], p[9], p[11], p[13], p[15]), lty = c(1,1,1,1,1,1,1,1))
  plot(PS.ROC15[,3], PS.ROC15[,2], "l", xlab="1-Specificity", ylab="Sensitivity", col=p[1], xlim=c(0,0.6), ylim=c(0,1), main = "")
  lines(MI.ROC15[,3], MI.ROC15[,2], col=p[3])
  points(CMA.ROC15[,3], CMA.ROC15[,2], col=p[6], pch=16)
  lines(hen.ROC15[,3], hen.ROC15[,2], col=p[7])
  points(log.ROC15[,3], log.ROC15[,2], col=p[9], pch=16)
  lines(RGS.ROC15[,3], RGS.ROC15[,2], col=p[11])
  lines(RS.ROC15[,3], RS.ROC15[,2], col=p[13])
  lines(HSMM.ROC15[2:19,3], HSMM.ROC15[2:19,2], col=p[15])
  legend("bottomright",c("PS", "MI", "CMA", "hen", "logisi", "RGS", "RS", "HSMM"),  bty="n", col=c(p[1], p[3], p[6], p[7], p[9], p[11], p[13], p[15]), lty = c(1,1,1,1,1,1,1,1))
  text(grconvertX(c(0.04, 0.51, 0.04, 0.51), from='ndc'),
       grconvertY(c(0.97, 0.97, 0.47, 0.47), from='ndc'),
       c('A',  'B',  'C',  'D'), xpd=NA, cex=1.5, font=2)
  dev.off()
  
  best.roc<-function(roc.res){
    c(roc.res[which.min(roc.res[,3]^2+(1-roc.res[,2])^2),], min(roc.res[,3]^2+(1-roc.res[,2])^2))
  }
  
  393 402  83.78600  10 5.41665 0.60185000  1
  
  s16$allb[[1]][38,]<-c(393, 406,   s15$spikes[[1]][393]-s15$spikes[[1]][392],  14,  s15$spikes[[1]][406]-s15$spikes[[1]][393],  (s15$spikes[[1]][406]-s15$spikes[[1]][393])/13,1) 
  
  res<-rbind(res, c(0,0,0))
  for (i in 12:14) {
    st<-s$spikes[[i]]
    start.true<- st[s$allb[[i]][,"beg"]]
    end.true<- st[s$allb[[i]][,"end"]]
    #si.adj<-subset(si.bursts[[i]], si.bursts[[i]][,"SI"]>j)
    si.adj<-mi.bursts[[i]]
    start.method<- st[si.adj[,"beg"]]
    end.method<- st[si.adj[,"end"]]
    res<-rbind(res, evaluate(st, start.true, end.true, start.method, end.method))
  }
  res<-rbind(res, c(0,0,0))
  for (i in 16:26) {
    st<-s$spikes[[i]]
    start.true<- st[s$allb[[i]][,"beg"]]
    end.true<- st[s$allb[[i]][,"end"]]
    #si.adj<-subset(si.bursts[[i]], si.bursts[[i]][,"SI"]>j)
    si.adj<-mi.bursts[[i]]
    start.method<- st[si.adj[,"beg"]]
    end.method<- st[si.adj[,"end"]]
    res<-rbind(res, evaluate(st, start.true, end.true, start.method, end.method))
  }
  
  
  ranksurprisemethod(s9$spikes[[1]], -log(0.05))
  
  ##
  # This function evaluates the performance of a burst detection method by comparing 
  # its output with the known location of bursts.
  #
  # Arguments:
  # spike.train = vector of spike timings
  # start.true = correct start times of bursts 
  # end.true = correct end times of bursts 
  # start.method = start times of bursts proposed by a burst detection algorithm
  # end.method = end times of bursts proposed by a burst detection algorithm
  #
  # Returns:
  # A list of the true positive rate, false positive rate and number of bursts detected
  ##
  
  
  evaluate <- function(spike.times, start.true, end.true, start.method, end.method){
    
    # Everything zero if no bursts were found by the method (marked with -1)
    if (max(start.method)<0){
      true.pos.rate <- 0
      false.pos.rate <- 0
      no.bursts <- 0
      
      # Otherwise evaluate proceeds
    }else{
      # Use spike.times and true start/end times to label spikes burst (=1) vs. no burst (=0)
      true.ID <- rep(0, length(spike.times))
      for (x in 1:length(start.true)){
        burst <- which(spike.times <= end.true[x] & spike.times >= start.true[x])
        true.ID[burst] <- 1
      }
      
      # Use spike.times and method start/end times to label spikes burst vs. no burst
      method.ID <- rep(0, length(spike.times))
      for (x in 1:length(start.method)){
        burst <- which(spike.times <= end.method[x] & spike.times >= start.method[x])
        method.ID[burst] <- 1
      }
      
      # Compare method.ID against true.ID
      total.pos <- length(which(true.ID==1))
      total.neg <- length(which(true.ID==0))
      score <-(2*true.ID)-method.ID       #true.pos=1, true.neg=0, false.pos=-1, false.neg=2
      true.pos <- length(which(score==1))
      false.pos <- length(which(score==-1))
      true.pos.rate <- true.pos/total.pos
      false.pos.rate <- false.pos/total.neg
      no.bursts <- length(start.method)
    }
    return(list(true.pos.rate, false.pos.rate, no.bursts))
  }
  
  ##Altered hennig method
  hennig.burst<-function(st) {
    bursts<-NULL
    allisi<-diff(st)
    isi.rank<-rank(allisi) #1 is smallest ISI
    st.length<-ceiling(max(st))
    spike.counts<-NULL
    for (i in 0:st.length-1) {
      spike.counts[i+1]<-sum((st>=i)*(st<(i+1)))
    }
    sc.hist<-hist(spike.counts, nclass=200)
    p.dist<-1-cumsum(sc.hist$counts/sum(sc.hist$counts))
    cutoff.indx<-sum(p.dist>0.05) ##Alter this
    theta.c<-max(c(2, ceiling(sc.hist$mids[cutoff.indx])))
    theta.c.off<-theta.c*0.5
    isi.rel.rank<-isi.rank/max(isi.rank)
    
    #t<-st[2:length(st)]#Why two??
    t<-st
    j<-1
    burst.on<-0
    bc<-1
    dt<-1
    burst.time<-NULL
    burst.end<-NULL
    burst.dur<-NULL
    burst.size<-NULL
    while (j<length(allisi)-theta.c) {
      if (burst.on==0 && isi.rel.rank[j]<0.5) {
        if (t[j+theta.c]<t[j]+dt){
          burst.on<-1
          burst.time[bc]<-t[j]
          brc<-j
        }
      } else if (burst.on==1) { 
        if (t[j+theta.c.off]>t[j]+dt) {
          burst.end[bc]<-t[j]
          burst.dur[bc]<-t[j]-burst.time[bc]
          burst.size[bc]<-j-brc
          bc<-bc+1
          burst.on<-0
        }
      }
      j<-j+1
    }
    if (burst.on==1) {
      tmp<-t[j]-burst.time[bc]
      burst.end[bc]<-burst.time[bc]+tmp
      burst.dur[bc]<-t[j]-burst.time[bc]
      burst.size[bc]<-j-brc
      bc<-bc+1
    }
    
    bursts<-cbind(burst.time, burst.end)
  }
  
  ##Altered PS method
  hennig.burst<-function(st) {
    bursts<-NULL
    allisi<-diff(st)
    isi.rank<-rank(allisi) #1 is smallest ISI
    st.length<-ceiling(max(st))
    spike.counts<-NULL
    for (i in 0:st.length-1) {
      spike.counts[i+1]<-sum((st>=i)*(st<(i+1)))
    }
    sc.hist<-hist(spike.counts, nclass=200)
    p.dist<-1-cumsum(sc.hist$counts/sum(sc.hist$counts))
    cutoff.indx<-sum(p.dist>0.05)
    theta.c<-max(c(2, ceiling(sc.hist$mids[cutoff.indx])))
    theta.c.off<-theta.c*0.5
    isi.rel.rank<-isi.rank/max(isi.rank)
    
    #t<-st[2:length(st)]#Why two??
    t<-st
    j<-1
    burst.on<-0
    bc<-1
    dt<-1
    burst.time<-NULL
    burst.end<-NULL
    burst.dur<-NULL
    burst.size<-NULL
    while (j<length(allisi)-theta.c) {
      if (burst.on==0 && isi.rel.rank[j]<0.5) {
        if (t[j+theta.c]<t[j]+dt){
          burst.on<-1
          burst.time[bc]<-t[j]
          brc<-j
        }
      } else if (burst.on==1) { 
        if (t[j+theta.c.off]>t[j]+dt) {
          burst.end[bc]<-t[j]
          burst.dur[bc]<-t[j]-burst.time[bc]
          burst.size[bc]<-j-brc
          bc<-bc+1
          burst.on<-0
        }
      }
      j<-j+1
    }
    if (burst.on==1) {
      tmp<-t[j]-burst.time[bc]
      burst.end[bc]<-burst.time[bc]+tmp
      burst.dur[bc]<-t[j]-burst.time[bc]
      burst.size[bc]<-j-brc
      bc<-bc+1
    }
    
    bursts<-cbind(burst.time, burst.end)
  }
  
  
  