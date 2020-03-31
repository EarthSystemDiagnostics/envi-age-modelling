# Process data for use in Bacon ----------

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


# Using only 14C dates at the moment. We need to harmonise the other dates
# and include them all. Some sites have a mix of date types


terr.14C.dat <- all.terr.dat %>%
  filter(age.type == "Radiocarbon years BP") %>%
  group_by(DataName) %>%
  # count the number of dates per site
  mutate(n.dates = n()) %>%
  ungroup()

length(unique(terr.14C.dat$DataName))


data.for.bacon.1 <- terr.14C.dat %>% 
  select(sample.id, DataName, depth, age, e.older) %>% 
  rename(sigma.age = e.older) %>%  # check whether this is OK
  # fix zero or missing sigma.age
  group_by(DataName) %>% 
  # first set NA values to zero
  mutate(sigma.age = ifelse(is.na(sigma.age), 0, sigma.age)) %>% 
  # then set zero values to mean of other sigma.age for that site
  # should probably take age into account here as uncertainty increases with age
  mutate(sigma.age = ifelse(sigma.age == 0, mean(sigma.age), sigma.age)) %>% 
  ungroup()


# Add a depth and date at surface

data.for.bacon <- data.for.bacon.1 %>%
  select(DataName) %>%
  distinct() %>%
  mutate(depth = 0, age = 0, 
         sigma.age = 100, # age uncertainty at fake surface points set to 100 years
         sample.id = 0,
         added.surface = TRUE) %>% 
  # inline bind
  bind_rows(data.for.bacon.1, .) %>%
  mutate(added.surface = ifelse(is.na(added.surface), FALSE, TRUE))

