rm(list=ls(all=TRUE))
set.seed(1)

nobs=1000
nseg=10
tmp=runif(nseg)
prob=tmp/sum(tmp); prob
partition=rmultinom(1,size=nobs,prob=prob)
mu=runif(nseg,min=5,max=50)

seg.index=rep(1:nseg,times=partition)
obs=data.frame(obs=rnorm(nobs,mean=mu[seg.index],sd=1),
               time1=1:nobs)

plot(obs~time1,data=obs)

setwd('U:\\GIT_models\\git_segmentation_model')
write.csv(obs,'fake data.csv',row.names=F)
