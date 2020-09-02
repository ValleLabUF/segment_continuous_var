gibbs.time.seg=function(dat,ngibbs,mu0,tau2,var) {  #where var is name of col of interest
  uni.id=unique(dat$id)
  dat.comp=dat[complete.cases(dat),]
  
  #to store results
  res.brks=vector("list", ngibbs)
  res.LML=matrix(NA,1,(ngibbs+1))
  res.nbrks=matrix(NA,1,(ngibbs+1))
  store.param=matrix(NA,ngibbs,2)
  
  #useful stuff
  min.time=min(dat.comp$time1)
  max.time=max(dat.comp$time1)
  
  #starting values
  breakpt=floor(mean(dat.comp$time1))
  
  
  for (i in 1:ngibbs){
    vals=samp.move(breakpt=breakpt,min.time=min.time,max.time=max.time,dat=dat.comp,tau2=tau2,
                   mu0=mu0,var=var)  
    breakpt=vals[[1]]
    
    #store results
    res.brks[[i]]<- breakpt
    store.param[i,]=c(length(vals[[1]]), vals[[2]])  # nbrks and LML
  }
  
  tmp=store.param[,1]
  res.nbrks[1,]=c(uni.id,tmp)
  colnames(res.nbrks)<- c('id', paste0("Iter_",1:ngibbs))
  
  tmp=store.param[,2]
  res.LML[1,]=c(uni.id,tmp)
  colnames(res.LML)<- c('id', paste0("Iter_",1:ngibbs))
  
  list(breakpt=res.brks, nbrks=res.nbrks, LML=res.LML)
}




#----------------------------------------------------
segment_time_continuous=function(data, ngibbs, mu0, tau2, var) {
  
  tic()  #start timer
  mod<- future_map(data, function(x) gibbs.time.seg(dat = x, ngibbs = ngibbs, mu0 = mu0,
                                                    tau2 = tau2, var = var), .progress = TRUE)
  toc()  #provide elapsed time
  
  
  brkpts<- map(mod, 1)  #create list of all sets breakpoints by ID
  
  nbrks<- purrr::map_dfr(mod, 2) %>%
    unlist() %>%
    matrix(.data, nrow = length(mod), ncol = (ngibbs + 1), byrow = T) %>%
    data.frame()  #create DF of number of breakpoints by ID
  names(nbrks)<- c('id', paste0("Iter_", 1:ngibbs))
  ncol.nbrks<- ncol(nbrks)
  nbrks<- nbrks %>%
    dplyr::mutate_at(2:ncol.nbrks, as.character) %>%
    dplyr::mutate_at(2:ncol.nbrks, as.numeric) %>%
    dplyr::mutate_at(1, as.character)
  
  LML<- purrr::map_dfr(mod, 3) %>%
    unlist() %>%
    matrix(.data, nrow = length(mod), ncol = (ngibbs + 1), byrow = T) %>%
    data.frame()  #create DF of LML by ID
  names(LML)<- c('id', paste0("Iter_", 1:ngibbs))
  ncol.LML<- ncol(LML)
  LML<- LML %>%
    dplyr::mutate_at(2:ncol.LML, as.character) %>%
    dplyr::mutate_at(2:ncol.LML, as.numeric) %>%
    dplyr::mutate_at(1, as.character)
  
  
  list(brkpts = brkpts, nbrks = nbrks, LML = LML)
}


