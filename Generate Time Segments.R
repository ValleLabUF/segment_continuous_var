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

# #smooth time series of R2n w/ moving average and add column for time
# dat.list<- map(dat.list, function(x) x %>% mutate(R2n.ma = as.numeric(forecast::ma(R2n, order = 10,
#                                                                         centre = TRUE)),
#                                                   time1 = 1:length(x)))

# #only select necessary cols
# dat.list<- map(dat.list, ~dplyr::select(., c(id, R2n.ma, time1)))


# Add time1 col
dat.list<- get_NSD(dat.list) %>% 
  map(., ~mutate(., time1 = 1:nrow(.)))

# Normalize NSD from 0 - 1
dat.list<- map(dat.list, ~mutate_at(., "NSD", ~{. / max(NSD)}))

# Add small value to keep obs where NSD = 0 from turning to -Inf after log-transform
dat.list<- map(dat.list, ~mutate_at(., "NSD", ~{. + 1e-09}))

#Log-transform NSD
dat.list<- map(dat.list, ~mutate(., log.NSD = log(.$NSD)))



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
                                  var = "log.NSD")
###Takes 2.5 min to run for 10000 iterations for all IDs


## Traceplots
#type is either 'nbrks' or 'LML' for y-axis label
identity<- names(dat.list)

traceplot(data = dat.res$nbrks, type = "nbrks", identity = identity)
traceplot(data = dat.res$LML, type = "LML", identity = identity)



## Determine maximum likelihood (ML) for selecting breakpoints
ML<- apply(dat.res$LML, 1, function(x) getML(dat = x, nburn = 5000))
brkpts<- getBreakpts(dat = dat.res$brkpts, ML = ML, identity = identity)


## Plots of data and breakpoints
plot.brks(data = dat.list, brkpts = brkpts, dat.res = dat.res, var = "log.NSD")


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
  geom_path(data = dat12 %>% filter(tseg == 1), aes(color=time1), size=0.5) +
  theme_bw() +
  scale_color_viridis_c() +
  coord_equal()




## Run a second pass of the model on log.NSD normalized by tseg

dat.list2<- map(dat.list, assign.time.seg, brkpts = brkpts) %>% 
  map(~bayesmove::df_to_list(.x, ind = "tseg")) %>%
  modify_depth(2, ~mutate(.x, log.NSD_n = log.NSD - min(log.NSD, na.rm = T))) %>%
  modify_depth(1, ~map_dfr(.x, `[`))



#priors
tau2=1
mu0=0

ngibbs = 10000

## Run Gibbs sampler
plan(multisession)  #run all MCMC chains in parallel
#select "multiprocess" if Unix or macOS & "multisession" if Windows
#refer to future::plan() for more details

dat.res2<- segment_time_continuous(data = dat.list2, ngibbs = ngibbs, mu0 = mu0, tau2 = tau2,
                                  var = "log.NSD_n")
###Takes 4 min to run for 10000 iterations for all IDs


## Traceplots
#type is either 'nbrks' or 'LML' for y-axis label
traceplot(data = dat.res2$nbrks, type = "nbrks", identity = identity)
traceplot(data = dat.res2$LML, type = "LML", identity = identity)


## Determine maximum likelihood (ML) for selecting breakpoints
ML2<- apply(dat.res2$LML, 1, function(x) getML(dat = x, nburn = 5000))
brkpts2<- getBreakpts(dat = dat.res2$brkpts, ML = ML2, identity = identity)


## Plots of data and breakpoints
plot.brks(data = dat.list2, brkpts = brkpts2, dat.res = dat.res2, var = "log.NSD_n")



dat_out2<- map(dat.list2, assign.time.seg, brkpts = brkpts2) %>% map_dfr(`[`)




### Map time segments of net-squared displacement

dat12<- dat_out2 %>% filter(id=="SNIK 12")

ggplot(data = dat12, aes(x,y)) +
  geom_path(aes(color=tseg), size=0.5, alpha=0.7) +
  theme_bw() +
  scale_color_viridis_c() +
  coord_equal()

ggplot() +
  geom_path(data = dat12, aes(x, y, color=time1), size=0.5) +
  theme_bw() +
  scale_color_viridis_c() +
  coord_equal() +
  facet_wrap(~tseg)






#create DF for breakpoints
#index brkpts for particular id
ind<- which(names(dat.list) == "SNIK 12")
breakpt<- brkpts2[ind,-1] %>%
  purrr::discard(is.na) %>%
  t() %>%
  data.frame()
names(breakpt)<- "breaks"


#plot comparison of NSD v log(NSD) in relation to breakpoints
ggplot() +
  geom_line(data = foo %>% filter(id == "SNIK 12"), aes(time1, NSD/max(NSD), color = "NSD")) +
  geom_line(data = foo %>% filter(id == "SNIK 12"), aes(time1, log(NSD)/max(log(NSD)),
                                                        color = "log(NSD)")) +
  geom_vline(data = breakpt, aes(xintercept = breaks)) + 
  facet_wrap(~id, scales = "free") +
  theme_bw() +
  labs(x = "Time", y = "Normalized Net Squared Displacement") +
  theme(strip.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 14))



#Export
write.csv(dat_out, "Snail Kite Gridded Data_TOHO_R2n.csv", row.names = F)













