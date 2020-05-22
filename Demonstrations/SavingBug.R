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

TestEm <- BuildNewEmulators(tData, HowManyEmulators = 2)
save(TestEm, file="TestEm.RData")
py_save_object(TestEm$mogp,"TestEm_mogp")
rm(TestEm)
load("TestEm.RData")
TestEm$mogp <- py_load_object("TestEm_mogp")
newDesign <- 2*randomLHS(100,3)-1
preds <- TestEm$mogp$predict(newDesign, deriv=FALSE)

####FIXED: New saving functions in BuildEmulator.R See implementation below####


TestEm <- BuildNewEmulators(tData, HowManyEmulators = 2)
save_ExUQmogp(TestEm, filename = "TestEm")
rm(TestEm)
TestEm <- load_ExUQmogp("TestEm")

newDesign <- 2*randomLHS(100,3)-1
preds <- TestEm$mogp$predict(newDesign, deriv=FALSE)

