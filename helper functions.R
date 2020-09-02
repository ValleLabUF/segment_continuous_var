assign.time.seg=function(dat,brkpts){  #add tseg assignment to each obs
  tmp=which(unique(dat$id) == brkpts$id)
  breakpt<- brkpts[tmp,-1] %>% purrr::discard(is.na) %>% as.numeric(.[1,])
  breakpt1=c(0,breakpt,Inf)
  n=length(breakpt1)
  res=matrix(NA,nrow(dat),1)
  for (i in 2:n){
    ind=which(breakpt1[i-1]<=dat$time1 & dat$time1<breakpt1[i])
    res[ind,]=i-1
  }
  dat$tseg<- as.vector(res)
  dat
}
#------------------------------------------------
df.to.list=function(dat) {  #only for id as col in dat
    id<- unique(dat$id)
    n=length(id)
    dat.list<- vector("list", n)
    names(dat.list)<- id
    
    for (i in 1:length(id)) {
      dat.list[[i]]<- dat[dat$id==id[i],]
    }
    dat.list
}
#------------------------------------------------
traceplot=function(data, type, identity) {  #create traceplots for nbrks or LML for all IDs
  for (i in 1:length(identity)) {
    par(ask=TRUE)
    plot(x=1:ngibbs, y=data[i,-1], type = "l", xlab = "Iteration",
         ylab = ifelse(type == "nbrks", "# of Breakpoints", "Log Marginal Likelihood"),
         main = paste("ID",identity[i]))
  }
  on.exit(par(ask = FALSE))
}
#---------------------------------------------
getML=function(dat,nburn) {  #select ML value that is beyond burn-in phase
  if (which.max(dat[-1]) < nburn) {
    ML<- dat[-1] %>% order(decreasing = T) %>% subset(. > nburn) %>% first()
  } else {
    ML<- which.max(dat[-1])
  }
  return(ML)
}
#---------------------------------------------
getBreakpts=function(dat,ML,identity) {  #extract breakpoints of ML per ID
  tmp<- list()
  
  for(i in 1:length(dat)) {
    ind<- ML[i]
    tmp[[i]]<- dat[[i]][[ind]]
  }
  
  names(tmp)<- identity
  max.length<- max(sapply(tmp, length))
  tmp<- lapply(tmp, function(x) { c(x, rep(NA, max.length-length(x)))})
  tmp<- map_dfr(tmp, `[`) %>% t() %>% data.frame()
  tmp<- cbind(id = identity, tmp)
  names(tmp)<- c('id', paste0("Brk_",1:(ncol(tmp)-1)))
  
  tmp
}
#------------------------------------------------
plot.brks.indiv=function(data, brkpts, var) {
  
  ind=which(unique(data$id) == brkpts$id)
  breakpt<- brkpts[ind,-1] %>% 
    discard(is.na) %>% t() %>% 
    data.frame()
  names(breakpt)<- "breaks"
  
  
  plot(data[,var], type="l", xlab = "Time", ylab = var, main = paste("ID",
                                                                 unique(data[,"id"])))
  abline(v=brkpts[ind,], col="grey")
  
}
#------------------------------------------------
plot.brks=function(dat.list, brkpts, var) {  
  
    par(ask = TRUE)
    map(dat.list, ~plot.brks.indiv(., brkpts = brkpts, var = var))
    par(ask = FALSE)
  
}
#----------------------------------
#Calc Net Squared Displacement
get_NSD_internal = function(dat) {
  
  # identify starting locs
  x0<- dat[1,"x"]
  y0<- dat[1,"y"]
  
  # calculate net squared displacement
  displ<- sqrt((dat[,"x"] - x0)^2 + (dat[,"y"] - y0)^2)
  nsd<- displ^2
  
  # add NSD to data
  dat$NSD<- nsd
  
  return(dat)
}

#----------------------------------
get_NSD = function(dat.list) {
  
  # map function to calculate NSD 
  tmp<- map(dat.list, get_NSD_internal)
  
  return(tmp)
}
