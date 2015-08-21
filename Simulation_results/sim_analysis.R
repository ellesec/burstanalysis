library(gridExtra)
library(parallel)

load("sim_data.RData")
source("~/sim_plot_functions.R")
set.seed(122)
types<-c("ps", "mi", "cma", "rs", "hennig", "rgs", "logisi.pasq", "clust", "hsmm")
cutoffs<-list(-log(0.01), NULL, FALSE, -log(0.01), NULL, NULL, NULL, 0.25, c(0.5, 0.9))
mi.par<-list(beg.isi=0.17,end.isi=0.3, min.ibi=0.2, min.durn=0.01, min.spikes=3)


sim.res$non.bursting<-sim.result(nb)
sim.res$non.stationary<-sim.result(ns)
sim.res$reg.bursting<-sim.result(sim.data$reg.bursting)
sim.res$long.bursts<-sim.result(sim.data$long.bursts)
sim.res$high.freq<-sim.result(sim.data$high.freq)
sim.res$noisy.bursts<-sim.tp.result(sim.data$noisy.bursts)


p<-plot.sib(sim.res[[1]], 0)
dummy <- p + theme(legend.position=c(0.2 ,0.2), legend.direction="vertical")
dummy.legend <- g_legend(dummy)
pdf("feat1.pdf", height=15, width=22)
grid.arrange(p, plot.sib(sim.res[[2]], 0), plot.sib(sim.res[[3]], 100), nrow=2, ncol=2, widths=c(5,5))
pushViewport(viewport(0.80, 0.4, 0.2, 0.5))
grid.draw(dummy.legend); popViewport()
pushViewport(viewport())
grid.text(LETTERS[1:3], gp=gpar(fontsize=20),
          x=unit(c(0.02, 0.52, 0.02), "npc"),
          y=unit(c(0.99,0.99, 0.5), "npc"))
dev.off()


pdf("feat2.pdf", height=21, width=22)
q<-plot.tp.fp(sim.res[[6]])
grid.arrange(plot.sib(sim.res[[4]], 100), plot.nbursts(sim.res[[4]], sim.data[[4]], 1),  plot.sib(sim.res[[5]], 100), plot.nbursts(sim.res[[5]], sim.data[[5]], 1), q[[1]], q[[2]], nrow=3, ncol=2)
grid.text(LETTERS[1:6], gp=gpar(fontsize=20),
          x=unit(c(0.02, 0.52, 0.02, 0.52, 0.02, 0.52), "npc"),
          y=unit(c(0.99,0.99, 0.66, 0.66, 0.33, 0.33), "npc"))
dev.off()

set.seed(123)
N<-round(runif(6, 1, 100))
pdf("sim_egs.pdf", height=15, width=14)
par(mfrow = c(6,1), mar=c(4,8,4,3), oma=c(0,0,0,0))
plot.sim.spikes(sim.data[[1]][[N[1]]]$spks, 231, 291)
lines(x=c(231, 236), y=c(-0.04, -0.04), lwd=3)
mtext(expression(bold("D5: Non-bursting   ")),2,4)
plot.sim.spikes(sim.data[[2]][[N[2]]]$spks, 240, 300)
mtext(expression(bold("D6: Non-stationary    ")),2,4)
plot.sim.spikes(sim.data[[3]][[N[3]]]$spks, 240, 300)
mtext(expression(bold("D7: Short bursts     ")),2,4)
plot.sim.spikes(sim.data[[4]][[N[4]]]$spks, 240, 300)
mtext(expression(bold("D8: Long bursts     ")),2,4)
plot.sim.spikes(sim.data[[5]][[N[5]]]$spks, 240, 300)
mtext(expression(bold("D9: High frequency       ")),2,4)
plot.sim.spikes(sim.data[[6]][[N[6]]]$spks, 240, 300)
mtext(expression(bold("D10: Noisy train      ")),2,4)
dev.off()