## Run for all cores with > 9 radiocarbon dates

library(tidyverse)
library(envibacon)

# load and filter data -----

all.terr.dat <- readr::read_csv2("inst/extdata/Dating_Data.csv") %>%
  filter(complete.cases(age, depth, e.older)) %>%
  mutate(sigma.age = e.older)

all.terr.14C.dat <- all.terr.dat %>%
  filter(age.type == "Radiocarbon years BP") %>%
  group_by(DataName) %>%
  mutate(n.dates = n()) %>%
  ungroup()

# Add a depth and date at surface
surface.ages <- all.terr.14C.dat %>%
  ungroup() %>%
  select(DataName) %>%
  distinct() %>%
  mutate(age = 0, depth = 0, sigma.age = 100, sample.id = 0,
         added.surface = TRUE)

all.terr.14C.dat.2 <- bind_rows(all.terr.14C.dat, surface.ages) %>%
  ungroup() %>%
  mutate(added.surface = ifelse(is.na(added.surface), FALSE, TRUE)) %>%
  #select(DataName, added.surface, everything()) %>%
  group_by(DataName) %>%
  mutate(n.dates = median(n.dates, na.rm = TRUE)) %>%
  ungroup()


terr.14C.min10.dates <- all.terr.14C.dat.2 %>%
  filter(n.dates >= 10)


## Get parameters for the runs

pars.df <- read.csv("inst/extdata/terr.14C.bacon.pars.csv", stringsAsFactors = FALSE)

pars.df <- pars.df %>%
  mutate(thick = round(100 / yrs_per_cm, 2),
         acc.mean = round(yrs_per_cm, 2),
         d.by = thick,
         K = (max.depth - min.depth) / thick) %>%
  rename(d.min = min.depth, d.max = max.depth) %>%
  tbl_df() %>%
  select(DataName, d.min, d.max, d.by, thick, acc.mean, K) %>%
  filter(DataName %in% terr.14C.min10.dates$DataName)

hist(pars.df$K)

# Some will take a long time.

MakeBaconDirs(terr.14C.min10.dates,
              filename = "terr_14C_min10_dates.csv",
              path = "inst/extdata/",
              suffix = "date.time",
              site.id = "DataName", sample.id = "sample.id",
              age = "age", age.err = "sigma.age", depth = "depth")


CreateParametersFiles("inst/extdata/terr_14C_min10_dates-2020.03.04_15-19-42/", pars.df)

# Transfer directory stucture and files to server
#
# scp -r "inst/extdata/terr_14C_min10_dates-2020.03.04_15-19-42/" adolman@linux4.awi-potsdam.de:rdata

#devtools::install_github("earthsystemdiagnostics/envi-age-modelling/envibacon")

#RunBaconDirs("rdata/terr_14C_min10_dates-2020.03.04_15-19-42/", runname = "")




