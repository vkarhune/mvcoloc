# rm(list=ls())

library(devtools)
library(roxygen2)

setwd("../")
create_package("mvcoloc")

use_mit_license()

use_r("simulate_phenotypes.R")
