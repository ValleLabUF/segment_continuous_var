#Analyze Snail Kite Data with Bayesian Partitioning Model

library(tidyverse)
library(progress)
library(furrr)
library(tictoc)
library(viridis)

source('gibbs functions.R')
source('helper functions.R')
source('gibbs sampler.R')

set.seed(1)

###############################
#### Load and Prepare Data ####
###############################

dat<- read.csv("Snail Kite Gridded Data_TOHO.csv", header = T, sep = ",")
dat$id<- as.character(dat$id)
dat.list<- bayesmove::df_to_list(dat = dat, ind = "id")


# Add time1 col
dat.list<- get_NSD(dat.list) %>% 
  map(., ~mutate(., time1 = 1:nrow(.)))

# Normalize NSD from 0 - 1
dat.list<- map(dat.list, ~mutate_at(., "NSD", ~{. / max(NSD)}))

# Add small value to keep obs where NSD = 0 from turning to -Inf after log-transform
dat.list<- map(dat.list, ~mutate_at(., "NSD", ~{. + 1e-12}))

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
###Takes 5 min to run for 10000 iterations for all IDs


## Traceplots
bayesmove::traceplot(data = dat.res$nbrks, ngibbs = ngibbs, type = "nbrks")
bayesmove::traceplot(data = dat.res$LML, ngibbs = ngibbs, type = "LML")



##Determine maximum a posteriori (MAP) estimate for selecting breakpoints
MAP.est<- bayesmove::get_MAP(dat.res$LML, nburn = ngibbs/2)
brkpts<- bayesmove::get_breakpts(dat = dat.res$brkpts, MAP.est = MAP.est)



## Plots of data and breakpoints
plot.brks(dat.list = dat.list, brkpts = brkpts, var = "log.NSD")






## Run a second pass of the model on log.NSD normalized by tseg

dat.list2<- map(dat.list, assign.time.seg, brkpts = brkpts) %>% 
  map(~bayesmove::df_to_list(.x, ind = "tseg")) %>%
  modify_depth(2, ~mutate(.x, log.NSD_n = log.NSD - min(log.NSD, na.rm = T))) #%>%
  # modify_depth(1, ~map_dfr(.x, `[`))

foo<- flatten(dat.list2)
names1<- map(dat.list2, length) %>% 
  unlist()
names(foo)<- paste(rep(names(names1), names1), names(foo), sep = "_")


#priors
tau2=1
mu0=0

ngibbs = 10000

## Run Gibbs sampler
plan(multisession)  #run all MCMC chains in parallel
#select "multiprocess" if Unix or macOS & "multisession" if Windows
#refer to future::plan() for more details

dat.res2<- segment_time_continuous(data = foo, ngibbs = ngibbs, mu0 = mu0, tau2 = tau2,
                                  var = "log.NSD_n")
###Takes 11 min to run for 10000 iterations for all segments


## Traceplots
#type is either 'nbrks' or 'LML' for y-axis label
bayesmove::traceplot(data = dat.res2$nbrks, ngibbs = ngibbs, type = "nbrks")
bayesmove::traceplot(data = dat.res2$LML, ngibbs = ngibbs, type = "LML")


##Determine maximum a posteriori (MAP) estimate for selecting breakpoints
MAP.est2<- bayesmove::get_MAP(dat.res2$LML, nburn = ngibbs/2)
brkpts2<- bayesmove::get_breakpts(dat = dat.res2$brkpts, MAP.est = MAP.est2)
brkpts2$id<- sub("\\_.*", "", brkpts2$id)  #modify id to remove segment index
brkpts2_merged<- brkpts2 %>%
  gather(key, value, -id) %>% 
  dplyr::select(-key) %>% 
  na.omit() %>% 
  pivot_wider(names_from = id, values_from = value, values_fn = list(value = list)) %>% 
  flatten() %>% 
  map(., . %>% 
        t() %>% 
        data.frame()) %>% 
  bind_rows() %>% 
  mutate(id = names(dat.list2)) %>% 
  # relocate(id)
  dplyr::select(id, everything())

## Plots of data and breakpoints
plot.brks(dat.list = dat.list, brkpts = brkpts2_merged, var = "NSD")





#create DF for breakpoints
#index brkpts for particular id
ind<- which(names(dat.list) == "SNIK 12")
breakpt<- brkpts[ind,-1] %>%
  purrr::discard(is.na) %>%
  t() %>%
  data.frame() %>% 
  rename(breaks = SNIK.12)

breakpt2<- brkpts2_merged[ind,-1] %>%
  purrr::discard(is.na) %>%
  t() %>%
  data.frame() %>% 
  rename(breaks = X22)

#plot comparison of NSD in relation to breakpoints
ggplot() +
  geom_line(data = dat.list$`SNIK 12`, aes(time1, NSD)) +
  geom_vline(data = breakpt, aes(xintercept = breaks), color = "darkturquoise") +
  # geom_vline(data = breakpt2, aes(xintercept = breaks), color = "firebrick") +
  facet_wrap(~id, scales = "free") +
  theme_bw() +
  labs(x = "Time", y = "NSD") +
  theme(strip.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 14))


######################################
#### Assign Spatial Time Segments ####
######################################

dat_out<- bayesmove::assign_tseg(dat.list, brkpts = brkpts)  #assign time seg and make as DF
dat_out$NSD<- dat$R2n  #include NSD on original scale

### Map time segments of net-squared displacement

dat12<- dat_out %>% filter(id=="SNIK 12")

ggplot(data = dat12, aes(x,y)) +
  geom_path(aes(color=tseg), size=0.5, alpha=0.7) +
  theme_bw() +
  scale_color_viridis_c() +
  coord_equal()


ggplot(data = dat12, aes(x,y)) +
  geom_path(size=0.75) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  coord_equal() +
  facet_wrap(~tseg)








#Export
# write.csv(dat_out, "Snail Kite Gridded Data_TOHO_R2n.csv", row.names = F)
# write.csv(brkpts, "Snail Kite NSD Breakpoints.csv", row.names = F)
