#MWE of mogp load problem:
#To be run with Working directory /ExeterUQ_MOGP/Demonstrations
mogp_dir <- "~/Dropbox/BayesExeter/mogp_emulator"
newwd <- getwd()
setwd("..")
source("BuildEmulator/BuildEmulator.R")
setwd(newwd)
load("ConvectionModelExample.Rdata")
TestEm <- BuildNewEmulators(tData, HowManyEmulators = 2)
print(TestEm[[1]])
save(TestEm, file="TestEm.RData")
rm(TestEm)
load("TestEm.RData")
print(TestEm[[1]])
