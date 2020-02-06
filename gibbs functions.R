get.summary.stats=function(breakpt,dat,var){
  breakpt1=c(0,breakpt,Inf)
  n=length(breakpt1)
  res=matrix(NA,n-1,3)
  colnames(res)=c('sum.x2','sum.x','n')
  for (i in 2:n){
    ind=which(breakpt1[i-1] <= dat$time1 & dat$time1 < breakpt1[i])
    tmp=dat[ind,]
    sum.x2=sum(tmp[,var]^2)
    sum.x=sum(tmp[,var])
    n=length(ind)
    res[i-1,]=c(sum.x2,sum.x,n)
  }
  res
}
#---------------------------------------------
log.marg.likel=function(tau2,mu0,summary.stats){
  #get sig2
  inv.sig2=summary.stats[,'n']+(1/tau2)
  sig2=1/inv.sig2
  
  #get mu1
  num1=summary.stats[,'sum.x']+(1/tau2)*mu0
  den1=summary.stats[,'n']+(1/tau2)
  mu1=num1/den1
  
  #get pieces of equation
  p1=summary.stats[,'n']*log(2*pi)
  p2=log(sig2/tau2)
  p3=(mu1^2)/sig2
  p4=summary.stats[,'sum.x2']+((mu0^2)/tau2)
  
  sum((1/2)*(-p1+p2+p3-p4))
}
#---------------------------------------------
samp.move=function(breakpt,min.time,max.time,dat,tau2,mu0,var){
  breakpt.old=breakpt
  p=length(breakpt)
  rand1=runif(1)	
  p0=1
  new.brk=sample((min.time+1):(max.time-1),size=1)  #don't propose a brkpt at beginning or end
  
  if (p == 1) {
    #birth
    if (rand1 < 1/2){
      breakpt.new=sort(unique(c(breakpt.old,new.brk)))
      p0=2/3 #death prob 2 -> 1 is (1/3) and birth prob 1 -> 2 is 1/2. 
    }
    #swap
    if (rand1 > 1/2) breakpt.new=new.brk
  }
  if (p > 1) {
    #birth
    if (rand1 < 1/3) {
      breakpt.new=sort(unique(c(breakpt.old,new.brk)))
    }
    #death
    if (rand1 > 1/3 & rand1 < 2/3) {
      ind=sample(1:length(breakpt.old),size=1)
      breakpt.new=breakpt.old[-ind]
      if (p==2) p0=3/2 #birth prob from 1 -> 2 is 1/2 and death prob from 2 -> 1 is 1/3
    }
    #swap
    if (rand1 > 2/3) {
      ind=sample(1:length(breakpt.old),size=1)
      breakpt.new=sort(unique(c(breakpt.old[-ind],new.brk)))
    }
  }
  
  #get sufficient statistics
  stats.old=get.summary.stats(breakpt=breakpt.old,dat=dat,var=var)
  stats.new=get.summary.stats(breakpt=breakpt.new,dat=dat,var=var)
  
  #get marginal loglikel
  pold=log.marg.likel(tau2=tau2,mu0=mu0,summary.stats=stats.old)
  pnew=log.marg.likel(tau2=tau2,mu0=mu0,summary.stats=stats.new)+log(p0)
  prob=exp(pnew-pold)
  rand2=runif(1)
  
  if (rand2<prob) return(list(breakpt.new, pnew))
  return(list(breakpt.old, pold))
}
