#Applies mi.find.bursts from sjemea to single spike train
MI.method<- function(spike.train){
  burst<-mi.find.bursts(spike.train)
  if (dim(burst)[1]<1) {
    burst<-NA
  }
  burst
}