# tidyverse packages for general data exploration
library(tidyverse)
library(lme4)

all.terr.dat <- readr::read_csv2("inst/extdata/terr_agemodel_data/Dating_Data.csv") %>%
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


# all.terr.c14.sed.rates.lm <- all.terr.c14.dat.2 %>%
#   #filter(added.surface == FALSE) %>%
#   group_by(DataName, n.dates) %>%
#   do({
#     tryCatch({
#       #lm1 <- MASS::rlm(age~depth, data = .)
#       lm1 <- lm(age~depth, data = .)
#       #broom::tidy(lm1)
#
#       data.frame(yrs_per_cm = coef(lm1)[2])
#     },
#     error = function(e) {data.frame(yrs_per_cm = NA)}
#     )
#   }) %>%
#   mutate(
#     cm_per_kyr = 1000 * 1/yrs_per_cm
#   )


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


pollen.depths.1 <- readxl::read_excel("/Users/andrewdolman/Dropbox/Work/AWI/Data/terrestrial-age-models/Neotoma_C14_with Min_Max_Pollen_2.xlsx")

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
  mutate(max.depth = max(c(max_pollen, max.14C.depth), na.rm = TRUE),
         min.depth = min(c(min_pollen, min.14C.depth), na.rm = TRUE),
         age.mod.length = max.depth - min.depth,
         thick = 100 / yrs_per_cm,
         K = age.mod.length / thick,
         K.older = ifelse(yrs.older > 0, yrs.older / thick, 0),
         K.younger = ifelse(yrs.younger > 0, yrs.younger / thick, 0),
         p.K.ex = (K.older + K.younger) / K)


pollen.depths.2 %>%
  #filter(n.dates > 9) %>%
  select(DataName, starts_with("K"), p.K.ex) %>%
  gather(var, val, -DataName) %>%
  #filter(val > 0) %>%
  ggplot(aes(x = val)) +
  geom_histogram() +
  facet_wrap(~var, scales = "free")

pollen.depths.2 %>%
  filter(n.dates > 9) %>%
  select(DataName, starts_with("K"), p.K.ex) %>%
  gather(var, val, -DataName) %>%
 # filter(val > 0) %>%
  ggplot(aes(x = val)) +
  geom_histogram() +
  facet_wrap(~var, scales = "free")

terr.14C.bacon.pars <- pollen.depths.2 %>%
  select(DataName, min.depth, max.depth, min.14C.depth, max.14C.depth, yrs_per_cm, n.dates) %>%
  mutate(yrs_per_cm = yrs_per_cm)

write.csv(terr.14C.bacon.pars, file = "inst/extdata/terr_agemodel_data/terr.14C.bacon.pars.csv",
          row.names = FALSE, quote = FALSE)

pollen.depths.2 %>%
  filter(n.dates > 9) %>%
  filter(K < 3*n.dates) %>%
  select(DataName, K, n.dates)

