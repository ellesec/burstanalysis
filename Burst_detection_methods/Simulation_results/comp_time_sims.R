set.seed(126)
period<-300
rate<-1
start<-0
pois.comp.time<-list()
for (i in 1:100){
pois.comp.time[[i]]<-sort(runif(rpois(1,period*rate), min=start, max=start+period))
}

types<-c("ps", "mi", "cma", "rs", "hennig", "rgs", "logisi", "hsmm")
cutoffs<-list(NULL, NULL, FALSE, -log(0.01), NULL, NULL, NULL, 0.5)
times<-system.time(sapply(pois.comp.time, function(x) burst.run(x, types[[1]], cutoffs[[1]])))
comp.time.df<-data.frame(val=times[3], type=types[[1]])
for (j in 2:7) {
  times<-system.time(sapply(pois.comp.time, function(x) burst.run(x, types[[j]], cutoffs[[j]])))
  comp.time.df<-rbind(comp.time.df, data.frame(val=times[3], type=types[[j]]))
}


qplot(factor(type), y=val, data=comp.time.df, geom="bar", stat="identity", fill=factor(type))
save(file="comp_time_sims.RData", pois.comp.time, comp.time.df)


burst.run<- function(spike.train, method, cutoff=NULL) {
  bursts <- switch(method, mi = MI.method(spike.train), 
                   ps = PS.method(spike.train), logisi = logisi.method(spike.train),
                   rs = RS.method(spike.train, cutoff), rgs=RGS.method(spike.train), hsmm=HSMM.method(spike.train, cutoff), 
                   cma=CMA.method(spike.train, brs.incl=cutoff), hennig=hennig.method(spike.train), stop(method, " : no such method for burst analysis"))
  bursts
}