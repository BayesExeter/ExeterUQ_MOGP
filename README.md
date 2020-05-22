# ExeterUQ_MOGP
An R interface to performing UQ with mogp_emulator.

This repository acts as an R front end for mogp_emulator (released in Python), which is available here: [https://github.com/alan-turing-institute/mogp_emulator](https://github.com/alan-turing-institute/mogp_emulator)

Specifically, we are using the implementation of mogp_emulator that allows for custom prior specification, custom mean function specification and for MAP estimation rather than maximum likelihood. **Currently this is only available on the devel branch of mogp_emulator**, so until this is updated **please use devel and not master.**

The goal of this repository is to enable users of UQ to easily and automatically fit emulators, produce diagnostics and history match, benefitting both from the fast GPU-optimised code offered by MOGP, and, using our defaults, from our experience in fitting emualtors to many outputs simultaneously, across a variety of contexts and in a safe way. Our interface is as simple to use as popular packages for emulation such as `Dicekriging`, and uses custom priors and MAP estimation to avoid common issues with maximum likelihood estimation (currently the default of mogp). Emulators, priors and mean functions are customisable so that, if you want to go your own way, ExeterUQ_MOGP can simply be the R interface to mogp_emulator that you need.

The pages on the website  [https://bayesexeter.github.io/ExeterUQ_MOGP/](https://bayesexeter.github.io/ExeterUQ_MOGP/)
offer a number of tutorials from using mogp_emulator directly (without our interface here), through to emulation and history matching at scale with our implementation. The Rmarkdown files for the tutorials are available in the demonstrations folder so that you can run the code for yourself.  
