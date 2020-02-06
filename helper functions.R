assign.time.seg=function(dat){  #add tseg assignment to each obs
  tmp=which(unique(dat$id) == brkpts$id)
  breakpt<- brkpts[tmp,-1] %>% discard(is.na) %>% as.numeric(.[1,])
  breakpt1=c(0,breakpt,Inf)
  tmp1<- which(diff(breakpt1) < 1)  #check for impossible time units
  breakpt1[tmp1+1]<- breakpt1[(tmp1+1)] + 1  #fix impossible time units
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
create.grid=function(dat,crs,extent,res,buffer) {
  grid<- raster(extent + buffer)
  res(grid)<- res
  proj4string(grid)<- crs
  grid[]<- 0
  
  grid
}
#------------------------------------------------
grid.summary.table=function(dat,crs,extent,res,buffer){  #dat must already have time.seg assigned
  
  #create grid and extract coords per cell
  grid<- create.grid(dat=dat, crs=crs, extent=extent, res=res, buffer=buffer)
  grid.cell.locs<- coordinates(grid) %>% data.frame()
  names(grid.cell.locs)<- c("x", "y")
  grid.cell.locs$grid.cell<- 1:length(grid)
  grid.coord<- grid.cell.locs[grid.cell.locs$grid.cell %in% dat$grid.cell,]
  
  
  grid.coord
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
         main = paste("ID",levels(identity)[data[i,1]]))
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
plot.brks.indiv=function(data, brkpts, dat.res, var) {
  
  ind=which(unique(data$id) == brkpts$id)
  breakpt<- brkpts[ind,-1] %>% discard(is.na) %>% t() %>% data.frame()
  names(breakpt)<- "breaks"
  
  
  plot(dat.list[[ind]][,var], type="l", ylab = var)
  abline(v=brkpts[ind,], col="grey")
  
  # print(
  #   ggplot(obs.long, aes(x=time, y=key, fill=value)) +
  #     geom_tile() +
  #     scale_fill_viridis_d("") +
  #     scale_y_continuous(expand = c(0,0)) +
  #     scale_x_continuous(expand = c(0,0)) +
  #     geom_vline(data = breakpt, aes(xintercept = breaks), color = viridis(n=9)[7], size = 0.35) +
  #     labs(x = "Observations", y = "Grid Cell") +
  #     theme_bw() +
  #     theme(axis.title = element_text(size = 18), axis.text = element_text(size = 16),
  #           title = element_text(size = 20)) +
  #     labs(title = paste("ID", unique(data$id)))
  # )
  
  
}
#------------------------------------------------
plot.brks=function(data, brkpts, dat.res, var) {  
  
    par(ask = TRUE)
    map(data, ~plot.brks.indiv(., brkpts = brkpts, dat.res = dat.res, var = var))
    par(ask = FALSE)
  
}
