## Process dating data to obtain parameters for Bacon runs.

# We need to:

# Get an estimate of the mean sediment accumulation rate for each site to use as the
# acc.mean parameter

# Get depth of the deepest point, including pollen samples, to know the depth range to model 


# Load data -----------

## Data should be in a directory called data-raw, inside the working directory


library(tidyverse)
library(broom)


all.terr.dat <- readr::read_csv2("data-raw/Dating_Data.csv") %>% 
  filter(complete.cases(age, depth)) %>%
  select(DataName, depth, age, e.older, e.young, age.type,
         material.dated, everything()) 

length(unique(all.terr.dat$DataName))

# 3145 unique sites


# Estimate acc.mean ----------

# Using only 14C dates at the moment. We need to harmonise the other dates
# and include them all. Some sites have a mix of date types


terr.14C.dat <- all.terr.dat %>%
  filter(age.type == "Radiocarbon years BP") %>%
  group_by(DataName) %>%
  # count the number of dates per site
  mutate(n.dates = n()) %>%
  ungroup()

length(unique(terr.14C.dat$DataName))


# Use a linear model to estimate acc.mean for each site

# Use the function do from dplyr to fit lm to each site and collect results
# in dataframe

acc.mean.ests <- terr.14C.dat %>% 
  group_by(DataName, n.dates) %>% 
  do({
    lm1 <- lm(age ~ depth, data = .)
    
    cfs <- broom::tidy(lm1)
    depth.cfs <- cfs[cfs$term == "depth", ]
    
    data.frame(
      yrs.per.cm = depth.cfs$estimate,
      SE.yrs.per.cm = depth.cfs$std.error
    )
  })

# check there is only one row per site
table(table(acc.mean.ests$DataName))


# take a look at the estimates

hist(acc.mean.ests$yrs.per.cm)

acc.mean.ests %>% 
  ggplot(aes(x = n.dates, y = yrs.per.cm)) +
  geom_point() +
  scale_x_continuous(trans = "sqrt", breaks = 2^(1:6))

# some extremely high, some negative, estimates are unreliable when there is 
# little data. This is why I previously used a hierarchical model, but that 
# introduced other problems.

# Obviously no estimate for locations with just one date

# can we use SE as estimate of reliability?

# "Hand" fit Bacon for those with low relative SE?
 


# Get deepest point ---------

# Need info about depths of pollen samples, using data file 
# "Neotoma_C14_with Min_Max_Pollen_2.xlsx" at the moment

pollen.depths.1 <- readxl::read_excel("data-raw/Neotoma_C14_with Min_Max_Pollen_2.xlsx")

# just need deepest (and maybe shallowest) point for each site

pollen.depths <- pollen.depths.1 %>%
  select(DataName, max_pollen, min_pollen) %>% 
  distinct()

# only 972 sites at this point

# Get the same from the dating data, use all including calendar but filter
# out lines with no age

dating.depths <- all.terr.dat %>% 
  select(DataName, depth, age) %>% 
  filter(complete.cases(age)) %>% 
  group_by(DataName) %>% 
  summarise(max.date.depth = max(depth),
            min.date.depth = min(depth))


# Join these together and get d.max and d.min 

depth.ranges <- full_join(pollen.depths, dating.depths, by = "DataName") %>% 
  group_by(DataName) %>% 
  mutate(d.max = max(max.date.depth, max_pollen, na.rm = T),
         d.min = min(min.date.depth, min_pollen, na.rm = T))

# we may well later set d.min to zero for all sites to model to the surface
 

## This is some analysis of how much we will be extrapolating

depth.range.analysis <- depth.ranges %>% 
  group_by(DataName) %>%
  mutate(length.dating = max.date.depth - min.date.depth,
         yrs.deeper = (max_pollen - max.date.depth),
         yrs.shallower = (min.date.depth - min_pollen),
         yrs.shallower.surface = (min.date.depth - 0),
         prop.deeper = yrs.deeper / length.dating,
         prop.shallower = yrs.shallower / length.dating,
         prop.shallower.surface = yrs.shallower.surface / length.dating)

# proportion of the age model that will be shallower than shallowest date
hist(depth.range.analysis$prop.shallower, 20)

# proportion of the age model that will be deeper than deepest date
hist(depth.range.analysis$prop.deeper, 20)


# Quite a lot of extrapolations going on.


# Create parameter dataframe to use with CreateParameterFiles() function -------


# join the acc.mean and d.min, d.max dataframes and select for required values

bacon.pars <- full_join(acc.mean.ests, depth.ranges, by = "DataName") %>% 
  rename(acc.mean = yrs.per.cm) %>% 
  select(DataName, acc.mean, d.max, d.min)


# add other parameters

# do not have info on PEAT or LAKE yet so keep mem parameter same for now
# But these parameters can be made conditional on type using ifelse statement
# or by joining another dataframe

bacon.pars <- bacon.pars %>% 
  mutate(acc.shape = 1.5,
         mem.mean = 0.7,
         mem.strength = 4,
         thick = 1, # set to 1 cm for now
         d.by = 1, # if d.by is 1, then the output will be at 1cm resolution
                  # even if the model (thick) is different
         d.min = 0 # overide d.min and set to 0 for all sites 
         )


### IMPORTANT #####

# Check the number of sections (K) that will be modelled.

bacon.pars.analysis <- bacon.pars %>% 
  mutate(K = (d.max - d.min) / thick)

bacon.pars.analysis %>% 
  ggplot(aes(x = K)) +
  geom_histogram()

# Some extremely large models!

# Should also check how many sections will be modelled between dating points

