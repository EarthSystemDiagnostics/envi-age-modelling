


#### Below here is code for running on server ########

## Change username and file paths accordingly.



## Bash code at terminal


# 1. Transfer the files you need to the server.

# In a terminal use scp to transfer files.

# You can either transfer the Dating_Data.csv, and also a file with parameters,
# e.g. terr.14C.bacon.pars.2.csv

# Or create a set of directories locally in R with MakeBaconDirs(...), 
# CreateParameterFiles(..). And transfer the whole top.level.dir structure


# scp -r "../working-data/terr_agemodel_data/Dating_Data.csv" adolman@linux4.awi-potsdam.de:rdata
# scp -r "../working-data/terr_agemodel_data/Some_Parameters.csv/" adolman@linux4.awi-potsdam.de:rdata

# OR e.g.

# scp -r "../working-data/terr_agemodel_data/terr_14C_min10_dates-2020.03.12_12-11-47/" adolman@linux4.awi-potsdam.de:rdata


# 2. In a new terminal log into the server with secure shell access

# ssh adolman@linux4.awi-potsdam.de 

# 3. In this terminal launch R by just typing "R"
# You should now be in an R Console


## R Code in Console in Terminal

# 4. Install envibacon

# devtools::install_github("earthsystemdiagnostics/envi-age-modelling/envibacon")
# you might need to install devtools first 
# install.packages("devtools")


# 5. Now you should be able to work in R on the server as if you are on your 
# local machine (but not RStudio, just the Console)

# RunBaconDirs(...) change the path etc.


### Bash code again in normal Terminal

# 6. Once it has finished you should be able to use the other terminal to 
# download the results to your machine with scp again.

# scp -r adolman@linux4.awi-potsdam.de:rdata/terr_14C_min10_dates-2020.03.12_12-11-47 "../working-data/terr_agemodel_data/terr_14C_min10_dates-2020.03.12_12-11-47_out/"


