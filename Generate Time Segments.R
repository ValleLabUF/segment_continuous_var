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
dat.list<- map(dat.list, function(x) x %>% mutate(R2n.ma = as.numeric(forecast::ma(R2n, order = 10,
                                                                        centre = TRUE)),
                                                  time1 = 1:length(x)))

# #only select necessary cols
# dat.list<- map(dat.list, ~dplyr::select(., c(id, R2n.ma, time1)))



#######################################
#### Run Gibbs Sampler for all IDs ####
#######################################

#priors
tau2=0.5
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
identity<- names(dat.list)

traceplot(data = dat.res$nbrks, type = "nbrks", identity = identity)
traceplot(data = dat.res$LML, type = "LML", identity = identity)



## Determine maximum likelihood (ML) for selecting breakpoints
ML<- apply(dat.res$LML, 1, function(x) getML(dat = x, nburn = 500))
brkpts<- getBreakpts(dat = dat.res$brkpts, ML = ML, identity = identity)


## Plots of data and breakpoints
plot.brks(data = dat.list, brkpts = brkpts, dat.res = dat.res, var = "R2n.ma")


######################################
#### Assign Spatial Time Segments ####
######################################



dat_out<- map(dat.list, assign.time.seg, brkpts = brkpts) %>% map_dfr(`[`)  #assign time seg and make as DF


### Map time segments of net-squared displacement

dat12<- dat_out %>% filter(id=="SNIK 12")

ggplot(data = dat12, aes(x,y)) +
  geom_path(aes(color=tseg), size=0.5, alpha=0.7) +
  theme_bw() +
  scale_color_viridis_c() +
  coord_equal()

ggplot(data = dat12, aes(x,y)) +
  geom_path(size=0.5, alpha=0.7, color = "grey75") +
  geom_point(data = dat12 %>% filter(tseg == 37), aes(color=time1), size=2) +
  theme_bw() +
  scale_color_viridis_c() +
  coord_equal()





#Export
write.csv(dat_out, "Snail Kite Gridded Data_TOHO_R2n.csv", row.names = F)

