# General outline of running Bacon in parallel


# 1. Process the Dating_Data and Pollen data (Neotoma_C14_with Min_Max_Pollen_2) 
# to get the parameters that will be used for each core. Output these parameters
# as a .csv file

# 2. Filter the Dating_Data for the sites you want to fit.

# 3. Use the function MakeBaconDirs() and CreateParameterFiles() to create a 
# set of folders, 1 for each site, containing a data file and parameters file

# 4. Transfer this set of folders to your user area on the server system using
# scp or whatever you know best

# 5. ssh into e.g. linux4@awi-potsdam.de

# 6. Launch an R session on the server

# 7. Install the package "envibacon" from Github if you haven't done this yet
# (or reinstall to update) devtools::install_github("")

# 8. Use function RunBaconDirs() to run Bacon for all the site subfolders

# 9. Transfer the output back to your machine with scp or similar

# 10. Use AggregateSummaryAgeModels() to collect all the age models together for 
# analysis

# 11. Stitch together all the pdfs into one file.


### Example code ------

if (!require("devtools")) {
  install.packages("devtools")
}

devtools::install_github("earthsystemdiagnostics/envi-age-modelling/envibacon")

library(tidyverse) 
library(envibacon)

# load and filter data -----------

# you will have to adjust paths for your system

all.terr.dat <- readr::read_csv2("../working-data/terr_agemodel_data/Dating_Data.csv") %>%
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
  # create a flag for the surface data point
  mutate(added.surface = ifelse(is.na(added.surface), FALSE, TRUE)) %>%
  group_by(DataName) %>%
  mutate(n.dates = median(n.dates, na.rm = TRUE)) %>%
  ungroup()


# Select which sites you want to run Bacon for -----------------

#terr.14C.min10.dates <- all.terr.14C.dat.2 %>%
#  filter(n.dates >= 10)


# Example with just a few sites

terr.14C.sample <- all.terr.14C.dat.2 %>%
  filter(DataName %in% c("DBATHTUB", "CAPPESJA"))


## Get parameters for this batch of Bacon runs -------------

### Parameters have been pre-calculated see get-bacon-pars.R

pars.df <- read.csv("../working-data/terr_agemodel_data/terr.14C.bacon.pars.2.csv", stringsAsFactors = FALSE)

# filter for just those sites we are going to use
pars.df <- pars.df %>%
  filter(DataName %in% terr.14C.sample$DataName)


# check distribution of number of sections
hist(pars.df$K)



# Make a set of folders for this run ------------

# This will create a top-level folder plus subfolders for each site

MakeBaconDirs(terr.14C.sample,
              filename = "terr.14C.sample.csv",
              path = "../working-data/terr_agemodel_data/",
              suffix = "date.time",
              site.id = "DataName", sample.id = "sample.id",
              age = "age", age.err = "sigma.age", depth = "depth")

# Add parameter files for each site

# You will need to change the path

CreateParameterFiles("../working-data/terr_agemodel_data/terr-2020.03.23_13-57-38/", pars.df)


# Run locally in Batch mode --------

## If you are on a Windows machine this will not actually run in parallel.

RunBaconDirs("../working-data/terr_agemodel_data/terr-2020.03.23_13-57-38/")


SummaryAgeMods <- AggregateSummaryAgeModels("../working-data/terr_agemodel_data/terr-2020.03.23_13-57-38/") %>% 
  tbl_df()

p <- SummaryAgeMods %>% 
  ggplot(aes(x = depth)) +
  geom_ribbon(aes(ymax = max, ymin = min), fill = "Grey") +
  geom_line(aes(y = mean)) +
  facet_wrap(~DataName)
p


## This function gets the age model specifically at the depths of the data by 
## interpolating between the modelled depths

AggregateAgeModelsAtDepths("../working-data/terr_agemodel_data/terr-2020.03.23_13-57-38/")




#### Below here is code for running on server ########


## Bash code at terminal

# Transfer directory stucture and files to server from a Terminal on your local machine

# scp -r "../working-data/terr_agemodel_data/terr_14C_min10_dates-2020.03.12_12-11-47/" adolman@linux4.awi-potsdam.de:rdata

# scp -r adolman@linux4.awi-potsdam.de:rdata/terr_14C_min10_dates-2020.03.12_12-11-47 "../working-data/terr_agemodel_data/terr_14C_min10_dates-2020.03.12_12-11-47_out/"



#devtools::install_github("earthsystemdiagnostics/envi-age-modelling/envibacon")

#RunBaconDirs("rdata/terr_14C_min10_dates-2020.03.12_12-11-47/", runname = "")




