
source("~/create_sims_functions.R")
#No bursts
set.seed(123)
period<-300
rate<-0.5
start<-0
end<-0
pois.no.bursts<-list()
for ( i in 1:50) {
  spike.times<-sort(runif(rpois(1,period*rate), min=start, max=start+period))
  isi<-diff(spike.times)
  rem.spikes<-which(isi<quantile(isi, 0.1))
  st<-NULL
  st$spks<-spike.times[-rem.spikes]
  pois.no.bursts[[i]]<-st
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
  temp<-NULL
  temp$spks<-gsum
  gamma.no.bursts[[i]]<-temp
}

all.no.bursts<-append(pois.no.bursts, gamma.no.bursts)



#Non stationary
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
  st<-NULL
  st$spks<-X2[-which(isi<q)]
  inhom.pois[[i]]<-st
}




#Regular bursts
set.seed(128)
rate<-0.2
period<-300
mean.length<-5

reg.bursts<-list()
for (i in 1:100){
  reg.bursts[[i]]<-pois.burst.noiseless(0.2, 300, 0)
}

#Long bursts
set.seed(129)
rate<-0.1
period<-300
mean.length<-18

long.bursts<-list()
for (i in 1:100){
  long.bursts[[i]]<-pois.burst.long(0.1, 300, 0)
}



#High frequency
set.seed(129)
rate<-1
period<-300
mean.length<-10

high.bursts<-list()
for (i in 1:100){
  high.bursts[[i]]<-pois.burst.high(1, 300, 0)
}




#Bursts with rate 0.5, length 8, on [-0.4, 0.4]. Noise gamma(0.5, 1)
set.seed(130)
rate<-0.5
period<-300
mean.length<-8

noisy.bursts<-list()
for (i in 1:100){
  noisy.bursts[[i]]<-pois.burst.noisy(0.5,300,0)
}


#Computational time
set.seed(126)
period<-300
rate<-1
start<-0
pois.comp.time<-list()
for (i in 1:100){
  pois.comp.time[[i]]<-sort(runif(rpois(1,period*rate), min=start, max=start+period))
}

sim.data<-list(all.no.bursts, inhom.pois, reg.bursts, long.bursts, high.bursts,  noisy.bursts, pois.comp.time)
names(sim.data)<-c("non.bursting", "non.stationary", "reg.bursting", "long.bursts", "high.freq", "noisy.bursts", "comp.time")