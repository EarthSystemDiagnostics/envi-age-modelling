# tidyverse packages for general data exploration
library(tidyverse)

# library(dplyr)
# library(tidyr)
# library(ggplot2)

library(lme4)

all.terr.dat <- readr::read_csv2("../working-data/terr_agemodel_data/Dating_Data.csv") %>%
  filter(complete.cases(age, depth, e.older)) %>%
  mutate(sigma.age = e.older)

all.terr.c14.dat <- all.terr.dat %>%
  filter(age.type == "Radiocarbon years BP") %>%
  group_by(DataName) %>%
  mutate(n.dates = n()) %>%
  ungroup()

# Add a depth and date at surface
surface.ages <- all.terr.c14.dat %>%
  ungroup() %>%
  select(DataName) %>%
  distinct() %>%
  mutate(age = 0, depth = 0, sigma.age = 100, sample.id = 0,
         added.surface = TRUE)

all.terr.c14.dat.2 <- bind_rows(all.terr.c14.dat, surface.ages) %>%
  ungroup() %>%
  mutate(added.surface = ifelse(is.na(added.surface), FALSE, TRUE)) %>%
  #select(DataName, added.surface, everything()) %>%
  group_by(DataName) %>%
  mutate(n.dates = median(n.dates, na.rm = TRUE)) %>%
  ungroup()



## Estimate mean sedimentation rate for each core to use as acc.mean in prior -----

# Use a hierarchical model to shrink extreme sed.rate estimates for cores with
# very few data points.

lmer1 <- all.terr.c14.dat.2 %>%
  filter(added.surface == FALSE) %>%
  lmer(age~ depth + (depth|DataName), data = .)

all.terr.c14.sed.rates <- coef(lmer1)$DataName %>%
  rownames_to_column(., var = "DataName") %>%
  as_tibble() %>%
  rename(yrs_per_cm = depth) %>%
  mutate(cm_per_kyr = 1000 * 1/yrs_per_cm) %>%
  select(-`(Intercept)`)

#
all.terr.c14.sed.rates <- all.terr.c14.dat.2 %>%
  select(DataName, n.dates) %>%
  distinct() %>%
  left_join(all.terr.c14.sed.rates, .) %>%
  mutate(n.date.cat =  cut(n.dates, c(0, 10, Inf)))


# Have a look at the sed.rates

all.terr.c14.sed.rates %>%
  #filter(cm_per_kyr > 0, cm_per_kyr < 50000) %>%
  ggplot(aes(x = cm_per_kyr)) +
  geom_histogram() +
  facet_wrap(~n.date.cat, scales = "free")


all.terr.c14.sed.rates %>%
  #filter(cm_per_kyr > 0, cm_per_kyr < 50000) %>%
  mutate(n.date.cat =  n.dates > 9) %>%
  ggplot(aes(x = yrs_per_cm)) +
  geom_histogram() +
  facet_wrap(~n.date.cat, scales = "free") +
  scale_x_continuous(trans = "sqrt", breaks = c(1,2,4,8,16, 32, 64, 128, 256))


## Some negative accumulation rates
# If negative set to bacon default of 20, which is actually close to overall mean.

all.terr.c14.sed.rates <- all.terr.c14.sed.rates %>%
  mutate(yrs_per_cm = ifelse(yrs_per_cm <= 0, 20, yrs_per_cm))

summary(all.terr.c14.sed.rates$yrs_per_cm)



## Get d.max and d.min from pollen data

n.dates <- all.terr.c14.dat %>%
  group_by(DataName) %>%
  summarise(n.dates = n(),
            max.14C.depth = max(depth, na.rm = TRUE),
            min.14C.depth = min(depth, na.rm = TRUE))


pollen.depths.1 <- readxl::read_excel("../working-data/terr_agemodel_data/Neotoma_C14_with Min_Max_Pollen_2.xlsx")

pollen.depths <- pollen.depths.1 %>%
  full_join(., n.dates) %>%
  filter(complete.cases(n.dates),
         n.dates > 1) %>%
  select(DataName, max_pollen, min_pollen, n.dates, max.14C.depth, min.14C.depth) %>%
  distinct() %>%
  left_join(all.terr.c14.sed.rates) %>%
  group_by(DataName) %>%
  mutate(length.dating = max.14C.depth - min.14C.depth,
         yrs.older = (max_pollen - max.14C.depth),
         yrs.younger = (min.14C.depth - min_pollen),
         perc.older = 100 * yrs.older / length.dating,
         perc.younger = 100 * yrs.younger / length.dating)

hist(pollen.depths$yrs.younger, 200)
hist(pollen.depths$yrs.older, 200)


pollen.depths.2 <-  pollen.depths %>%
  group_by(DataName) %>%
  mutate(min.14C.depth = 0) %>%
  mutate(max.depth = ceiling(max(c(max_pollen, max.14C.depth), na.rm = TRUE)),
         min.depth = floor(min(c(min_pollen, min.14C.depth), na.rm = TRUE)),
         age.mod.length = max.depth - min.depth,
         thick = 100 / yrs_per_cm,
         K = age.mod.length / thick)


# Put some sensible limits on K and thick and K.per.n.dates
terr.14C.bacon.pars <- pollen.depths.2 %>%
  mutate(d.by = thick,
         acc.mean = round(yrs_per_cm, 2)) %>%
  rename(d.min = min.depth, d.max = max.depth) %>%
  # mutate(yrs.per.slice = thick * acc.mean)
  mutate(K = ifelse(K > 200, 200, K),
         thick = d.max / K,
         K.per.n = K / n.dates) %>%
  mutate(thick = ifelse(thick > 20, 20, thick),
         K = d.max / thick,
         K.per.n = K / n.dates) %>%
  mutate(K.per.n = ifelse(K.per.n < 3, 3, K.per.n),
         K = K.per.n * n.dates,
         thick = d.max / K) %>%
  mutate(thick = ifelse(thick > 20, 20, thick),
         K = d.max / thick,
         K.per.n = K / n.dates,
         yrs.per.slice = thick * acc.mean)

hist(subset(terr.14C.bacon.pars, terr.14C.bacon.pars$n.dates > 8)$K)


terr.14C.bacon.pars.2 <- terr.14C.bacon.pars %>%
  mutate(thick = round(thick, 2),
         acc.mean = round(yrs_per_cm, 2),
         d.by = thick,
         K = round(K)) %>%
  select(DataName, d.min, d.max, d.by, thick, acc.mean, n.dates, K) %>%
  mutate(acc.shape = 1.5, mem.strength = 4, mem.mean = 0.7)


write.csv(terr.14C.bacon.pars.2, file = "../working-data/terr_agemodel_data/terr.14C.bacon.pars.2.csv",
          row.names = FALSE, quote = FALSE)

