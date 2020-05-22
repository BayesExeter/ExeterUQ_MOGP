# ExeterUQ_MOGP
An R interface to performing UQ with mogp_emulator

This repository acts as a front end for mogp, which is available here: [https://github.com/alan-turing-institute/mogp_emulator](https://github.com/alan-turing-institute/mogp_emulator)

Specifically, we are using the implementation that allows for custom prior specification, custom mean function specification and for MAP estimation rather than maximum likelihood. **Currently this is only available on the devel branch of mogp**, so until this is updated **please use devel and not master**

This repository offers our own custom default subjective prior specification, designed to be as similar as possible to the ExeterUQ code: [https://github.com/BayesExeter/ExeterUQ](https://github.com/BayesExeter/ExeterUQ) that is used by different groups for climate model tuning.

We provide methods for easy customisable emulation, diagnostics and history matching. Work through the demonstrations: first to see how mogp works through the ```mogp_emulator_demo.Rmd```, then build an ExeterUQ_mogp through ```Testing-ExeterUQ-mogp.Rmd```. Learn how diagnostics work with ```mogp_emulator_diagnostics.Rmd``` and see history matching through ```HistoryMatchingwithExeterUQ.Rmd```.

Go to [Vignettes](Vignettes)