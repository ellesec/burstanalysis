MI.method.array<-function(s){
  bursts<-spikes.to.bursts(s,"mi")
  null.burst<-which(sapply(bursts, function(x) dim(x)[1])<1)
  if (length(null.burst)){
    bursts[null.burst]<-NA
  }
  bursts
}

MI.method<- function(spike.train){
  burst<-mi.find.bursts(spike.train)
  if (dim(burst)[1]<1) {
    burst<-NA
  }
  burst
}