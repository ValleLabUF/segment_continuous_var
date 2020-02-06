#Analyze Snail Kite Data with Bayesian Partitioning Model

library(tidyverse)
library(progress)
library(furrr)
library(tictoc)
library(viridis)
library(forecast)

source('gibbs functions.R')
source('helper functions.R')
source('gibbs sampler.R')


###############################
#### Load and Prepare Data ####
###############################

dat<- read.csv("Snail Kite Gridded Data_TOHO.csv", header = T, sep = ",")
dat.list<- df.to.list(dat = dat)

#smooth time series of R2n w/ moving average and add column for time
dat.list<- map(dat.list, function(x) x %>% mutate(R2n.ma = ma(R2n, order = 10, centre = TRUE),
                                                  time1 = 1:length(x)))

#only select necessary cols
dat.list<- map(dat.list, ~dplyr::select(., c(id, R2n.ma, time1)))



#######################################
#### Run Gibbs Sampler for all IDs ####
#######################################

#priors
tau2=1
mu0=0

ngibbs = 10000

## Run Gibbs sampler
plan(multisession)  #run all MCMC chains in parallel
                    #select "multiprocess" if Unix or macOS & "multisession" if Windows
                    #refer to future::plan() for more details

dat.res<- segment_time_continuous(data = dat.list, ngibbs = ngibbs, mu0 = mu0, tau2 = tau2,
                                  var = "R2n.ma")
###Takes 8 min to run for 10000 iterations for all IDs


## Traceplots
#type is either 'nbrks' or 'LML' for y-axis label
identity<- unique(dat$id)

traceplot(data = dat.res$nbrks, type = "nbrks", identity = identity)
traceplot(data = dat.res$LML, type = "LML", identity = identity)



## Determine maximum likelihood (ML) for selecting breakpoints
ML<- apply(dat.res$LML, 1, function(x) getML(dat = x, nburn = 500))
brkpts<- getBreakpts(dat = dat.res$brkpts, ML = ML, identity = identity)


## Heatmaps
plot.brks(data = dat.list, brkpts = brkpts, dat.res = dat.res, var = "R2n.ma")


######################################
#### Assign Spatial Time Segments ####
######################################



dat_out<- map(dat.list, assign.time.seg) %>% map_dfr(`[`)  #assign time seg and make as DF
dat_out<- rbind(dat_out, dat.ex)  #bring back in excluded data occupying < 3 cells
dat_out<- dat_out[order(dat_out$id, dat_out$date),]  #reorder DF by id and date

setwd("~/Documents/Snail Kite Project/Data/R Scripts/activcenter_subset_locations")
write.csv(dat_out, "Snail Kite Gridded Data_TOHO.csv", row.names = F)

