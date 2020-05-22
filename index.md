# ExeterUQ_MOGP
An R interface to performing UQ with mogp_emulator

This repository acts as a front end for mogp, which is available here: [https://github.com/alan-turing-institute/mogp_emulator](https://github.com/alan-turing-institute/mogp_emulator)

Specifically, we are using the implementation that allows for custom prior specification, custom mean function specification and for MAP estimation rather than maximum likelihood. **Currently this is only available on the devel branch of mogp**, so until this is updated **please use devel and not master**

This repository offers our own custom default subjective prior specification, designed to be as similar as possible to the ExeterUQ code: [https://github.com/BayesExeter/ExeterUQ](https://github.com/BayesExeter/ExeterUQ) that is used by different groups for climate model tuning.

We provide methods for easy customisable emulation, diagnostics and history matching. The pages below offer a number of tutorials from using mogp_emulator directly (without our interface here), through to emulation and history matching at scale with our implementation. The Rmarkdown files for the tutorials are available in the demonstrations folder so that you can run the code for yourself.

- [Running mopg\_emulator by hand in R](mogp_emulator_demo)
- [Easy Fitting of Emulators with ExeterUQ_mogp](Testing_ExeterUQ_mogp.md)
- [Saving an ExeterUQ_mogp Emulator](SavingAndLoading)
- [Explaining the ExeterUQ_mogp Prior](Subjective-Prior-distributions)
- [Diagnostics with ExeterUQ_mogp](mogp_emulator_diagnostics)
- [History matching with ExeterUQ_mogp](HistoryMatchingwithExeterUQ-MOGP)