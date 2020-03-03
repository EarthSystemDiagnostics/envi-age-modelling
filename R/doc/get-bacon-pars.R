# tidyverse packages for general data exploration
library(tidyverse)

all.terr.dat <- readr::read_csv2("/Users/andrewdolman/Dropbox/Work/AWI/Data/terrestrial-age-models/Dating_Data.csv") %>%
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
  mutate(age = 0, depth = 0, sigma.age = 100,
         added.surface = TRUE)

all.terr.c14.dat.2 <- bind_rows(all.terr.c14.dat, surface.ages) %>%
  mutate(added.surface = ifelse(isTRUE(added.surface), TRUE, FALSE)) %>%
  group_by(DataName) %>%
  mutate(n.dates = median(n.dates, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(depth.c = depth - mean(depth),
         age.c = age - mean(age))



all.terr.c14.sed.rates <- all.terr.c14.dat.2 %>%
  #filter(added.surface == FALSE) %>%
  group_by(DataName, n.dates) %>%
  do({
    tryCatch({
      #lm1 <- MASS::rlm(age~depth, data = .)
      lm1 <- lm(age~depth, data = .)
      #broom::tidy(lm1)

      data.frame(yrs_per_cm = coef(lm1)[2])
    },
    error = function(e) {data.frame(yrs_per_cm = NA)}
    )
  }) %>%
  mutate(
    cm_per_kyr = 1000 * 1/yrs_per_cm
  )


## Estimate mean sedimentation rate for each core to use as acc.mean in prior

library(lme4)
lmer1 <- all.terr.c14.dat.2 %>%
  filter(added.surface == FALSE) %>%
  lmer(age~ depth + (depth|DataName), data = .)
lmer2 <- lmer(age~ 1  + (depth|DataName), data = all.terr.c14.dat.2)
summary(lmer1)
summary(lmer2)

hist(coef(lmer1)$DataName[,"depth"], breaks = 20)
hist(coef(lmer2)$DataName[,"depth"])

summary(lmer1)

anova(lmer1, lmer2)


all.terr.c14.sed.rates <- coef(lmer1)$DataName %>%
  rownames_to_column(., var = "DataName") %>%
  as_tibble() %>%
  rename(yrs_per_cm = depth) %>%
  mutate(cm_per_kyr = 1000 * 1/yrs_per_cm) %>%
  select(-`(Intercept)`)

all.terr.c14.sed.rates %>%
  #filter(cm_per_kyr > 0, cm_per_kyr < 50000) %>%
  ggplot(aes(x = cm_per_kyr)) +
  geom_histogram()

all.terr.c14.sed.rates %>%
  #filter(cm_per_kyr > 0, cm_per_kyr < 50000) %>%
  ggplot(aes(x = yrs_per_cm)) +
  geom_histogram()

all.terr.c14.sed.rates %>%
  left_join(all.terr.c14.dat.2, .) %>%
  #filter(cm_per_kyr > -40000, cm_per_kyr < 50000) %>%
  ggplot(aes(x = n.dates, y = yrs_per_cm)) +
  geom_point()


summary(all.terr.c14.sed.rates$yrs_per_cm)


## Get d.max and d.min

n.dates <- all.terr.c14.dat %>%
  group_by(DataName) %>%
  summarise(n.dates = n(),
            max.14C.depth = max(depth, na.rm = TRUE),
            min.14C.depth = min(depth, na.rm = TRUE),
            geo.mean.d.depth = (mean(sqrt(diff(sort(depth)))))^2)

all.terr.c14.dat %>%

  filter(DataName == "WESTHAWK",
         depth < 200) %>%
  ggplot(aes(x = depth, y = age, colour = lab.no, group = NA)) +
  geom_pointrange(aes(ymin = age - e.young, ymax = age + e.older)) +
  geom_smooth(method = "lm") +
  expand_limits(x = 0, y = 0)

pollen.depths.1 <- readxl::read_excel("/Users/andrewdolman/Dropbox/Work/AWI/Data/terrestrial-age-models/Neotoma_C14_with Min_Max_Pollen_2.xlsx")

pollen.depths <- pollen.depths.1 %>%
  left_join(., n.dates) %>%
  filter(complete.cases(n.dates),
         n.dates > 1) %>%
  select(DataName, max_pollen, min_pollen, n.dates, max.14C.depth, min.14C.depth, geo.mean.d.depth) %>%
  distinct() %>%
  left_join(all.terr.c14.sed.rates)


pollen.depths %>%
  #filter(max_pollen < 2e4) %>%
  ggplot(aes(x = max_pollen)) +
  geom_histogram()

pollen.depths.2 <-  pollen.depths %>%
  group_by(DataName) %>%
  mutate(max.depth = max(c(max_pollen, max.14C.depth)),
         min.depth = min(c(min_pollen, min.14C.depth)),
         age.mod.length = max.depth - min.depth,
         thick = 100 / yrs_per_cm,
         thick2 = geo.mean.d.depth / 3,
         K1 = age.mod.length / thick,
         K2 = age.mod.length / thick2,
         K3 = n.dates * 3,
         K = max(c(K1, K2, K3), na.rm = TRUE))

pollen.depths.2 %>%
  filter(K < 1000,
         n.dates > 9) %>%
  ggplot(aes(x = K, y = K2)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1)

pollen.depths.2 %>%
  #filter(K < 1000,
  #       n.dates > 9) %>%
  ggplot(aes(x = K, y = K3)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1)

pollen.depths.2 %>%
  #filter(K < 1000) %>%
  ggplot(aes(x = K)) +
  geom_histogram()

summary(pollen.depths.2$K)
hist(pollen.depths.2$age.mod.length / pollen.depths.2$K)

good.terr.c14.dat <- all.terr.c14.dat %>%
  left_join(., sed.rates) %>%
  filter(n.dates >= 10,
         cm_per_kyr > 0,
         cm_per_kyr < 1000)

length(unique(good.terr.c14.dat$DataName))

sort(unique(good.terr.c14.dat$DataName))

good.terr.c14.dat %>%
  ggplot(aes(x = depth, y = age, group = DataName)) +
  geom_line()


n.slices <- good.terr.c14.dat %>%
  group_by(DataName) %>%
  summarise(median.d.depth = median(diff(sort(depth))[diff(sort(depth))>0]),
            thick.med = median.d.depth / 1,
            core.length = diff(range(depth)),
            n.slices.med = core.length / thick.med,

            mean.d.depth = mean(diff(sort(depth))),
            sqrt.mean.d.depth = (mean(sqrt(diff(sort(depth)))))^2,

            # smallest non.zero diff
            geo.mean.d.depth = exp(mean(log(diff(sort(depth)) + min(diff(sort(depth))[diff(sort(depth)) > 0])))),

            thick.mean = mean.d.depth / 1,
            geo.thick.mean = geo.mean.d.depth / 1,
            n.slices.mean = core.length / thick.mean,
            n.slices.geo.mean = core.length / geo.thick.mean)

summary(n.slices$core.length)

n.slices %>%
  ggplot(aes(x = n.slices.mean, y = n.slices.geo.mean)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1)

n.slices %>%
  ggplot(aes(x = n.slices.geo.mean, y = n.slices.med)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1)

n.slices %>%
  ggplot(aes(x = core.length, y = n.slices.mean)) +
  geom_point()


n.slices %>%
  ggplot(aes(x = mean.d.depth, y = geo.mean.d.depth)) +
  geom_point()

n.slices %>%
  ggplot(aes(x = sqrt.mean.d.depth, y = geo.mean.d.depth)) +
  geom_point()+
  geom_abline(intercept = 0, slope = 1)



hist(n.slices$geo.thick.mean)
n.slices$thick.med
n.slices$median.d.depth


