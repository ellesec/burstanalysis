plot.bursts<-function(s, cell, beg, end) {
  spike.train<-s$spikes[[cell]]
  plot( c(beg, end), c(0,1), type='n', bty='n', yaxt="n", xaxt="n", xlab="", main="", ylab="")
  b<-s$allb[[cell]]
  ymi<-plot.jitter(dim(b)[1], 0.4)
  n<-length(spike.train)
  ys<-numeric(n)
  segments(spike.train, ys, spike.train, ys+0.8, lwd=0.2) 
  segments(spike.train[b[,"beg"]], ymi, spike.train[b[,"end"]], ymi,col="#336699", lwd=3)
}

pdf("ex_data.pdf", height=3, width=12)
par(mfrow = c(2,1), mar=c(1,1,1,1), oma=c(0,0,0,0))
plot.bursts(demas$s15, 12, 782.5, 842.5)
lines(x=c(782.5, 788), y=c(-0.04, -0.04), lwd=2)
plot.bursts(demas$s15, 22, 75,135)
dev.off()