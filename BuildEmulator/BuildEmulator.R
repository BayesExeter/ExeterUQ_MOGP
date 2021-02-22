packages <- c('reticulate','pracma', 'invgamma', 'GenSA', 'far', 'fields', 'lhs', 'maps', 'mco', 'mvtnorm', 'ncdf4', 'parallel', 'shape', 'tensor', 'withr', 'loo','MASS')
sapply(packages, require, character.only = TRUE, quietly = TRUE)
source("BuildEmulator/AutoLMcode.R")
source("BuildEmulator/CustomPredict.R")



twd <- getwd()
setwd(mogp_dir)
mogp_emulator <- import("mogp_emulator")
mogp_priors <- import("mogp_emulator.Priors")
mogp_kernels <- import("mogp_emulator.Kernel")
setwd(twd)

extract_targets <- function(df, target_cols = c("y")) {
#'@Description Separate a data frame into inputs, targets, and inputdict for use with MOGP class.
#'@param df A data frame with column names for the inputs and targets.
#'@param target_cols The names of the columns of the data frame that are targets.
#'@return A list with a matrix of inputs a matrix of targets and an inputdict that can be used with the MOGP class.
  for (t in target_cols) {
    stopifnot(t %in% names(df))#Make sure the targets are in the input data frame
  }
  n_targets <- length(target_cols)
  inputs <- matrix(NA, ncol=ncol(df) - n_targets, nrow=nrow(df))
  targets <- matrix(NA, ncol=nrow(df), nrow=n_targets)
  inputdict <- dict()
  
  input_count <- 1
  target_count <- 1
  
  for (n in names(df)) {
    if (n %in% target_cols) {
      targets[target_count,] <- as.matrix(df[n])
      target_count <- target_count + 1
    } else {
      inputs[,input_count] <- as.matrix(df[n])
      inputdict[n] <- as.integer(input_count - 1)
      input_count <- input_count + 1
    }
  }
  
  if (n_targets == 1) {
    targets <- c(targets)
  }
  return(list(inputs, targets, inputdict))
}

FormulaToString <- function(tformula){
  #'@description Parse and R formula into a string
  #'@param tformula an R formula of the form y~I(x1)+I(x_2)+ etc
  #'@return A string version to be passed into MOGP
  f = as.character(tformula)
  paste(f[2], f[1], f[3], sep="")
}

alphaFun <- function(x,M,V){
  #'@description A polynomial in the shape of an inverse gamma distribution with coefficients depending on the mode M and Variance V the roots of which are needed to compute the shape.
  #'@param x the shape of the inverse gamma (to be rooted via Newton Raphson)
  #'@param M the mode of the inverse gamma distribution
  #'@param V the variance of the inverse gamma distribution.
  V*x^3 - (4*V+M^2)*x^2 + (5*V-2*M^2)*x - (2*V+M^2)
}

alphaFun2 <- function(M, V){
  #'@param M the mode of the inverse gamma distribution
  #'@param V the variance of the inverse gamma distribution
  #'@description finding the real root > 2 that will be the shape of our inverse gamma prior. 
  #'@detail polyroot() finds the complex roots. Either there are 3 real roots (we choose the highest >2 so that a variance exists) or 1 real and 2 complex conjugate roots. In the latter case the real root will still have a tiny numerically 0 imaginary part, and we look to isolate this first.
  troots <- polyroot(c(-(2*V+M^2),5*V-2*M^2,-(4*V+M^2),V))
  treal <- which(round(Im(troots), digits=6) ==0)
  if(length(treal)>1){
    troots <- troots[treal]
    return(troots[which(Re(troots)>2)[1]])
  }
  else if(length(treal)<1)
    stop("No real roots for the shape parameter, check roots and consider rounding by fewer digits")
  else{
    myroot <- Re(troots[treal])
    if(!(myroot>2))
      stop("Only real root for shape equation is < 2 and prior has no variance")
    else
      return(myroot)
  }
}

invgamMode <- function(bound, tmode){
  #'@description Finds the shape and rate of an inverse gamma distribution given it's mode and a bound on the uncertainty.
  #'@param bound The largest value you want. The variance of the inverse gamma is set so that crossing this bound is a 6 standard deviation event
  #'@param tmode the mode of the inverse gamma you want to fit.
  #'@details The alpha parameter is found using newtonRaphson on the alphaFun function. Beta is calculated directly given alpha.
  #'@return A list with the alpha and beta parameters
  nstd <- bound-tmode
  varSig <- (nstd/6)^2 #6 can be changed to make bound crossing more or less likely.
  ###
  #New part finding alpha using polyroot
  #alphaSig <- newtonRaphson(alphaFun, x0=tmode, M=tmode, V=varSig,maxiter=1000)$root
  alphaSig <- alphaFun2(M=tmode,V=varSig)
  ###
  betaSig <- tmode*(alphaSig+1)
  #Uncomment the below 3 lines to see the distribution and the probability of exceeding the bound.
  #hist(rinvgamma(10000,alphaSig,betaSig),breaks=100)
  #%abline(v=tmode)
  #print(1-pinvgamma(bound,alphaSig,betaSig))
  list(alpha=alphaSig,beta=betaSig)
}

GetStarts <- function(lm.emulator, d, Choices){
  MeanStart <- as.vector(summary(lm.emulator$linModel)$coef[,1])
  tSig <- summary(lm.emulator$linModel)$sigma^2
  c(MeanStart, rep(0,d), sqrt((1-Choices$NuggetProportion)*tSig), sqrt(Choices$NuggetProportion*tSig))
}

GetPriors <- function(lm.emulator, d, Choices, ActiveVariables){
  #'@description This function constructs a subjective prior for the parameters of a GP emulator to be fit by MO_GP. 
  #'@param lm.emulator A lm emulator list (see AutoLMCode). The only required element of this list is the linModel component. lm.emulator$linModel is an lm object fitted to the target data. Custom lm objects can be specified using lm.emulator = list(linModel=lm(...), entering your own formulae). It is suggested that this is done elsewhere and passed here.
  #'@param d the number of input parameters
  #'@param Choices A list containing hyperprior choices that control the subjective prior.
  #'@param NonInformativeRegression If TRUE, a uniform prior is used for the regression parameters.
  #'@param NonInformativeCorrelationLengths If TRUE, a uniform prior is used for the correlation parameters. 
  #'@param NonInformativeSigma If TRUE, a uniform prior is used for the sigma parameter.
  #'@param NonInformativeNugget If TRUE, a uniform prior is used for the regression parameters.
  #'@param BetaRegressMean Prior mean for the regression coefficients. The intercept is given a uniform prior, all other regression terms get a Normal prior with mean BetaRegressMean and variance "BetaRegressSigma"
  #'@param BetaRegressSigma Prior variance for the regression coefficients.
  #'@param DeltaActiveMean The mean of a lognormal prior for the active inputs (see details). 
  #'@param DeltaActiveSigma The variance of a lognormal prior for the active inputs (see details).
  #'@param DeltaInactiveMean The mean of a lognormal prior for the inactive inputs (see details).
  #'@param DeltaInactiveSigma The variance of a lognormal prior for the inactive inputs (see details).
  #'@param NuggetProportion What proportion of the data variability is suspected to be nugget. Only used if Nugget="fit"
  #'@param Nugget either a number or "fit", "fixed" or "adaptive". This is seen by MOGP which will only fit the nugget (and hence require a prior) if "fit" is chosen. If the others, the nugget is fixed to the value given or fitted to avoid numerical instabilities adaptively. (See MOGP documentation.)
  #'@param ActiveVariables Indices indicating which parameters in the data are active
  #'@return A list of priors in the following order. A prior for the intercept (defaults to NULL), p priors for the regression coefficients, d priors for the correlation lengths, a prior for sigma squared and a prior for the nugget (NULL if Choices$Nugget != "fit")
  #'@details The linear mean has p coefficients not associated with the intercept. Each coefficient is Beta ~ Normal(BetaRegressMean, BetaRegressSigma). Our default is based on mean 0 variance 10 and is weakly informative.
  #'@details Input parameters are classified either as active or inactive. These definitions depend on whether the parameters were highlighted by preliminary fitting algorithms in our linear modelling code (see AutoLMcode) and whether the user asks for them specifically to be included (so note they may do nothing, but still get an "active" prior. There are 2 classes of prior we use for the correlation lengths. The first is for the "active" ones: log(delta) ~ Normal(DeltaActiveMean, DeltaActiveSigma) (so that delta is lognormal). default this to N(0,0.5) on the basis that we expect there to be correlation and so that they are not too long (we penalise the rigde on the likelihood surface). Inactive parameters get the same type of prior with very strong default values N(5,0.005), giving the whole distribtion at contributions to the correlation near 1. As our covariance functions are separable, the inactive variables therefore don't alter the fit. All values for these parameters are controllable by the user.
  #'@details. Sigma^2 ~ InvGamma(M,V) where we use a parameterisation of the inverse gamma based on the mode and variance. The mode is chosen so reflect the variance we can explain with simple linear fits and the variance is set using invgamMode with a bound based on the total variability in the original data. The idea behind the prior is that we allow Sigma^2 to be as large as the variance of the data so that our emulator explains none of the variability, but we expect to explain as much as could be explained by a preliminary fit of a basic model.
  #'@details The nugget distribution, if required, is found as with sigma but using the choice NuggetProportion to reflect what percentage of variability we might expect to be nugget. We may add the ability to add a user prior here. Only really important for internal variability models like climate models. Most deterministic models should use "adaptive" or "fixed" nugget. 
  p <- length(lm.emulator$linModel$coefficients)-1
  Priors <- lapply(1:(1+p+d+1+1), function(e) NULL)
  if(!is.null(Choices$intercept))
    print("NULL intercept fitted by default: change code")
  #Regression priors
  if(!(Choices$NonInformativeRegression)){
    Betas <- lapply(1:p, function(e) mogp_priors$NormalPrior(Choices$BetaRegressMean, Choices$BetaRegressSigma))
    for(i in 2:(p+1)){#the first element is NULL for the intercept term
    Priors[[i]] <- Betas[[i-1]]
    }
  }
  #Correlation length priors
  if(!(Choices$NonInformativeCorrelationLengths)){
    Deltas <- lapply(1:d, function(k) {if(k %in% ActiveVariables) 
      {mogp_priors$NormalPrior(Choices$DeltaActiveMean,Choices$DeltaActiveSigma)} 
      else {mogp_priors$NormalPrior(Choices$DeltaInactiveMean,Choices$DeltaInactiveSigma)}})
    for(j in (p+2):(p+1+d)){
      Priors[[j]] <- Deltas[[j-p-1]]
    }
  }
  #Sigma and nugget priors
  ModeSig <- var(lm.emulator$linModel$residuals)
  boundSig <- ModeSig/(1-summary(lm.emulator$linModel)$r.squared)
  if(!(Choices$NonInformativeSigma)){
    SigmaParams <- invgamMode((1-Choices$NuggetProportion)*boundSig,ModeSig)
    Sigma <- mogp_priors$InvGammaPrior(SigmaParams$alpha,SigmaParams$beta)
    Priors[[d+p+2]] <- Sigma
  }
  if(Choices$Nugget=="fit"){
    if(!Choices$NonInformativeNugget){
      NuggetParams <- invgamMode(Choices$NuggetProportion*boundSig ,Choices$NuggetProportion*ModeSig)
      Nugget <- mogp_priors$InvGammaPrior(NuggetParams$alpha, NuggetParams$beta)
      Priors[[d+p+3]] <- Nugget
    }
  }
    return(Priors)
}

GetKernel <- function(kernel.type = "Gaussian") {
  #'@description This function specifies a kernel for a GP emulator to be fit by MO_GP. 
  #'@param kernel.type A single sting corresponding to either Gaussian, or Materm kernels.
  #'@details A function returns a MO_GP object.
  if(kernel.type == "Gaussian") {
    Kernel <- mogp_kernels$SquaredExponential()
  }
  if(kernel.type == "Matern52") {
    Kernel <- mogp_kernels$Matern52()
  }
  return(Kernel)
}

#Below our default choices list for fitting subjective Bayesian emulators. Most terms are discussed in the function GetPriors. The last 3 control the linear model fitting.
# lm.tryFouriers is a Boolean that indicates whether our automatic linear model fitting code should explore fourier terms.
# lm.maxOrder specifies what degree of polynomials can be fitted in our automatic linear models code. 
# lm.maxdf specifies how many degrees of freedom should be spent on fitting a linear model. The default is ~10% and is used if set to NULL
choices.default <- list(NonInformativeRegression=FALSE, 
                        NonInformativeCorrelationLengths = FALSE,
                        NonInformativeSigma = FALSE,
                        NonInformativeNugget = FALSE,
                        DeltaActiveMean = -0.25, DeltaActiveSigma = 0.14,
                        DeltaInactiveMean=-5, DeltaInactiveSigma=0.005,
                        BetaRegressMean = 0, BetaRegressSigma = 100,
                        NuggetProportion=0.05, Nugget = "fit", 
                        lm.tryFouriers=FALSE, lm.maxOrder=NULL, 
                        lm.maxdf=NULL)

DegreesFreedomDefault <- function(ChoiceList, N){
#' @param N is number of data points. 
#' @param ChoiceList is a list of lists with lm.maxdf in each list
#' @description If NULL is specified the degrees of freedom to be used by the lm is calculated as about 10% of N + 1
#' @return a vector of degrees of freedom
  t.dfs <- rep(NA, length(ChoiceList))
  for(i in 1:length(ChoiceList)){
    if(is.null(ChoiceList[[i]]$lm.maxdf))
      t.dfs[i] <- ceiling(N/10)+1
    else
      t.dfs[i] <- ChoiceList[[i]]$lm.maxdf
  }
  return(t.dfs)
}
EMULATE.lm <- function(Response, tData, dString="tData",maxdf=NULL,tcands=cands,tcanfacs=canfacs,TryFouriers=FALSE,maxOrder=NULL){
  #' Generate a linear model for Gaussian Process (GP) Emulator
  #' 
  #' @param Response a string corresponding to the response variable
  #' @param tData a data frame containing the inputs and outputs
  #' @param dString a string corresponding to the name of a data frame
  #' @param maxdf a maximim number of degrees of freedom that we are expecting to 
  #' lose by adding terms to the linear model
  #' @param tcands a vector of parameter names
  #' @param tcanfacs a vector of parameter names that are factors
  #' @param TryFouriers a logical argument with default FALSE. TRUE allows
  #' the Fourier transformation of the parameters
  #' @param maxOrder maximum order of Fourier terms (Fourier series)
  #' 
  #' @return A linear model and the parameters allowing refitting (see AutoLMcode.R)
  #' 
  if(is.null(maxdf)){
    maxdf <- ceiling(length(tData[,1])/10)
  }
  startLM <- eval(parse(text=paste("lm(",Response,"~1,data=",dString,")",sep="")))
  if(maxdf < 2){
    return(list(linModel=startLM,Names=NULL,mainEffects=NULL,Interactions=NULL,Factors=NULL,FactorInteractions=NULL,ThreeWayInters=NULL,DataString=dString,ResponseString=Response,tData=tData,BestFourier=TRUE))
  }
  if(TryFouriers){
    msl <- list(linModel=startLM,Names=NULL,mainEffects=NULL,Interactions=NULL,Factors=NULL,FactorInteractions=NULL,ThreeWayInters=NULL,DataString=dString,ResponseString=Response,tData=tData,BestFourier=TRUE,maxOrder=maxOrder)
  }
  else{
    msl <- list(linModel=startLM,Names=NULL,mainEffects=NULL,Interactions=NULL,Factors=NULL,FactorInteractions=NULL,ThreeWayInters=NULL,DataString=dString,ResponseString=Response,tData=tData,BestFourier=FALSE)
  }
  added <- AddBest(tcands,tcanfacs,msl)
  for(i in 1:30){
    added <- AddBest(tcands,tcanfacs,added)
    if(!is.null(added$Break))
      break
  }
  print(summary(added$linModel))
  if(nrow(tData)-1-added$linModel$df > maxdf){
    NRM <- nrow(tData) - 1 - added$linModel$df.residual - maxdf
    trm <- removeNterms(N=NRM,linModel=added$linModel,dataString=added$DataString,responseString=added$ResponseString,Tolerance=NULL,Names=added$Names,mainEffects=added$mainEffects,Interactions=added$Interactions,Factors=added$Factors,FactorInteractions=added$FactorInteractions,ThreeWayInters=added$ThreeWayInters,tData=added$tData,Fouriers = added$Fouriers)
    print(summary(trm$linModel))
    trm$pre.Lists <- get.predict.mats(trm$linModel)
    trm$DataString <- dString
    trm$ResponseString <- Response
    NRM <- nrow(tData) - 1 - trm$linModel$df.residual - maxdf
    while(NRM>0){
      trm <- removeNterms(N=NRM, linModel=trm$linModel, dataString=trm$DataString, responseString=trm$ResponseString, Tolerance=NULL, Names=trm$Names,mainEffects=trm$mainEffects, Interactions=trm$Interactions, Factors=trm$Factors, FactorInteractions=trm$FactorInteractions, ThreeWayInters=trm$ThreeWayInters, tData=added$tData, Fouriers = trm$Fouriers)
      print(summary(trm$linModel))
      trm$pre.Lists <- get.predict.mats(trm$linModel)
      trm$DataString <- dString
      trm$ResponseString <- Response
      NRM <- nrow(tData) - 1 - trm$linModel$df.residual - maxdf
    }
    return(trm)
  }
  else{
    return(added)
  }
}

EMULATE.lm.old <- function(Response, tData, dString="tData",maxdf=NULL,tcands=cands,tcanfacs=canfacs,TryFouriers=FALSE,maxOrder=NULL){
  #' Generate a linear model for Gaussian Process (GP) Emulator
  #' 
  #' @param Response a string corresponding to the response variable
  #' @param tData a data frame containing the inputs and outputs
  #' @param dString a string corresponding to the name of a data frame
  #' @param maxdf a maximim number of degrees of freedom that we are expecting to 
  #' lose by adding terms to the linear model
  #' @param tcands a vector of parameter names
  #' @param tcanfacs a vector of parameter names that are factors
  #' @param TryFouriers a logical argument with default FALSE. TRUE allows
  #' the Fourier transformation of the parameters
  #' @param maxOrder maximum order of Fourier terms (Fourier series)
  #' 
  #' @return A linear model and the parameters allowing refitting (see AutoLMcode.R)
  #' 
  if(is.null(maxdf)){
    maxdf <- ceiling(length(tData[,1])/10)
  }
  startLM <- eval(parse(text=paste("lm(",Response,"~1,data=",dString,")",sep="")))
  if(maxdf < 2){
    return(list(linModel=startLM,Names=NULL,mainEffects=NULL,Interactions=NULL,Factors=NULL,FactorInteractions=NULL,ThreeWayInters=NULL,DataString=dString,ResponseString=Response,tData=tData,BestFourier=TRUE))
  }
  if(TryFouriers){
    msl <- list(linModel=startLM,Names=NULL,mainEffects=NULL,Interactions=NULL,Factors=NULL,FactorInteractions=NULL,ThreeWayInters=NULL,DataString=dString,ResponseString=Response,tData=tData,BestFourier=TRUE,maxOrder=maxOrder)
  }
  else{
    msl <- list(linModel=startLM,Names=NULL,mainEffects=NULL,Interactions=NULL,Factors=NULL,FactorInteractions=NULL,ThreeWayInters=NULL,DataString=dString,ResponseString=Response,tData=tData,BestFourier=FALSE)
  }
  added <- AddBest(tcands,tcanfacs,msl)
  for(i in 1:30){
    added <- AddBest(tcands,tcanfacs,added)
    if(!is.null(added$Break))
      break
  }
  print(summary(added$linModel))
  if(nrow(tData)-1-added$linModel$df > maxdf){
    trm <- removeNterms(N=500,linModel=added$linModel,dataString=added$DataString,responseString=added$ResponseString,Tolerance=sum(sort(anova(added$linModel)$"Sum Sq"))*(1e-4),Names=added$Names,mainEffects=added$mainEffects,Interactions=added$Interactions,Factors=added$Factors,FactorInteractions=added$FactorInteractions,ThreeWayInters=added$ThreeWayInters,tData=added$tData,Fouriers = added$Fouriers)
    print(summary(trm$linModel))
    trm2 <- removeNterms(N=max(c(0,length(tData[,1])-maxdf-1-trm$linModel$df)),linModel=trm$linModel,dataString=added$DataString,responseString=added$ResponseString,Tolerance=sum(sort(anova(added$linModel)$"Sum Sq"))*(5e-1),Names=trm$Names,mainEffects=trm$mainEffects,Interactions=trm$Interactions,Factors=trm$Factors,FactorInteractions=trm$FactorInteractions,ThreeWayInters=trm$ThreeWayInters,tData=added$tData,Fouriers=trm$Fouriers)
    print(summary(trm2$linModel))
    trm2$pre.Lists <- get.predict.mats(trm2$linModel)
    trm2$DataString <- dString
    trm2$ResponseString <- Response
    return(trm2)
  }
  else{
    return(added)
  }
}

BuildNewEmulators <- function(tData, HowManyEmulators, 
                              additionalVariables=NULL, 
                              Choices = lapply(1:HowManyEmulators,
                                               function(k) choices.default), meanFun = "linear", 
                              kernel = c("Gaussian"),...){
  #'@description Builds MO_GP emulators for a data frame 
  #'@param tData This is the data frame you wish to use. The format should be D+1+Q where Q is the number of targets each occupying one of the last Q columns, D is the number of inputs occupying the first D columns. The D+1th column should be a vector of random normal numbers (between -1 and 1) called "Noise". We use it in our LM code to stop our algorithms overfitting by selecting signals by chance. All variables should have names.
  #'@param HowManyEmulators How many emulators are required. The code will fit the first HowManyEmulators GPs up to Q.
  #'@param additionalVariables Are there variables that must be "active" (e.g. you really want to use them in a decision problem or similar) or should be included? Often all variables are specified here, but it defaults to NULL.
  #'@param Choices A list of choices with each of the HowManyEmulators elements being a choice list compatible with GetPriors() (see the documentation for GetPriors)
  #'@param meanFun Currently a single string either "fitted", or "linear" ("constant" to come). If "fitted", our custom global mean functions are fitted and then used to fit the GP. Recommended for higher dimensions and for history matching. If "linear", a linear mean function is fitted to all emulators and using additionalVariables. A list implementation will be considered in future versions. Could also make a list where it could be a formula (is.formula)
  #'@param kernel A vector of strings that corresponds to the type of kernel either "Gaussian", or "Matern52"
  #'Default is to use Gaussian kernel for all emulators.
  #'@details If mean functions are not given (an option that will be added soon) our automatic LM code fits a global mean function for each metric. Get Priors is then used to extract the priors before we establish an MOGP and fit the parameters by MAP estimation. MAP improves on MLE here as we avoid the ridge on the likelihood surface for GPs.
  #'@return A list with 2 elements. 1 the mogp, 2 a list containing the elements used for fitting: the mean functions (containing element linModel as the lm object), the design, a list of active input indices for each emulator (to be used for subsetting the design to produce diagnostics), and the prior choices.
  ###Mean function selection###
  lastCand <- which(names(tData)=="Noise")
  if(length(lastCand)<1)
    stop("tData should have a column called 'Noise' separating the inputs and outputs.")
  if(is.null(HowManyEmulators))
    HowManyEmulators <- length(names(tData)) - lastCand
  if(!(HowManyEmulators == length(names(tData)) - lastCand)){
    tData <- tData[,c(1:lastCand,(lastCand+1):(lastCand+HowManyEmulators))]
  }
  if(meanFun =="fitted"){
    tdfs <- DegreesFreedomDefault(Choices, N=length(tData[,1]))
    lm.list = lapply(1:HowManyEmulators, function(k) 
      try(EMULATE.lm(Response=names(tData)[lastCand+k],
                     tData=tData, tcands=names(tData)[1:lastCand],
                     tcanfacs=NULL, 
                     TryFouriers=Choices[[k]]$lm.tryFouriers, 
                     maxOrder=Choices[[k]]$lm.maxOrder,
                     maxdf = tdfs[k])))
  
    }
  else if(meanFun == "linear"){
    if(is.null(additionalVariables))
      stop("When specifying linear meanFun, please pass the active inputs into additionalVariables")
    linPredictor <- paste(additionalVariables,collapse="+")
    lm.list = lapply(1:HowManyEmulators, function(k) list(linModel=eval(parse(text=paste("lm(", paste(names(tData[lastCand+k]), linPredictor, sep="~"), ", data=tData)", sep="")))))
  }
  else{
    stop("meanFun must either be 'fitted' or 'linear' in this version")
  }
  ###Prepare the data for MOGP### 
  tfirst <- lastCand + 1
  target_names <- names(tData)[tfirst:length(names(tData))]
  target_list <- extract_targets(tData[,-which(names(tData)=="Noise")], target_names)
  inputs <- target_list[[1]]
  targets <- target_list[[2]]
  inputdict <- target_list[[3]]
  d <- dim(inputs)[2]
  
  if(meanFun=="fitted"){
    ActiveVariableIndices <- lapply(lm.list, function(tlm) which((names(tData)%in%additionalVariables) | (names(tData)%in%tlm$Names) | (names(tData) %in% names(tlm$Fouriers))))
  }
  else if(meanFun == "linear"){
    ActiveVariableIndices <- lapply(lm.list, function(tlm) which(names(tData)%in%additionalVariables))
  }
  ###Prepare the mean functions for MOGP### 
  mean_func.list.MGP <- lapply(lm.list, function(e) FormulaToString(formula(e$linModel)))
  
  ###Establish the priors for the emulators###
  Priors <- lapply(1:HowManyEmulators, function(l) GetPriors(lm.list[[l]], d=d, Choices[[l]], ActiveVariableIndices[[l]]))
  
  ###Establish the kernel types for MOGP###
  if(length(kernel) == 1) {
    Kernels <- lapply(1:HowManyEmulators, function(l) GetKernel(kernel))
  } else {
    Kernels <- lapply(1:HowManyEmulators, function(l) GetKernel(kernel[l])) 
  }
  ### Where to Start the optimisation ###
    Starts <- lapply(1:HowManyEmulators, function(l) GetStarts(lm.list[[l]], d=d, Choices[[l]]))
  
  ###Establish and fit the MOGP###
  Emulators <- mogp_emulator$MultiOutputGP(inputs, targets, mean = mean_func.list.MGP,
                                           priors = Priors, inputdict = inputdict,
                                           nugget = lapply(Choices,function(e) e$Nugget), 
                                           kernel = Kernels)
  Emulators <- mogp_emulator$fit_GP_MAP(Emulators, n_tries=1, ftol=1e-06, theta0 = Starts, maxiter=1000)
  
  ###Prepare return objects###
  Design <- tData[,1:(lastCand-1), drop = FALSE]
  fitting <- list(lm.object = lm.list,
                  Design = Design, 
                  ActiveIndices = ActiveVariableIndices,
                  PriorChoices = Choices)
  
  return(list(mogp = Emulators, # call mogp
              fitting.elements= fitting))
}

save_ExUQmogp <- function(ExUQmogp, filename){
  #'@param ExUQmogp An ExeterUQ_mogp object. This is a list with 2 elements: $mogp (an mogp object) and $fitting.elements (various objects used in the fitting and useful for prediction.)
  #'@param filename This is a string without an extension. Do not use extensions such as .RData
  #'@description the mogp part of the ExUQmogp is a python object, but ExUQmogp is an Rlist. This code saves the mogp in 2 parts: An RData file with name "filename.RData" and a python object with name "filename_mogp". These 2 files are then recombined with load_ExUQmogp(filename).
  if(!all(names(ExUQmogp) == c("mogp", "fitting.elements")))
    stop("Object to save must be an Exeter UQ mogp object: a list with 2 elements: the mogp and the fitting elements.")
  save(ExUQmogp, file = paste(filename,"RData",sep="."))
  py_save_object(ExUQmogp$mogp, paste(filename,"mogp",sep="_"))
}

load_ExUQmogp <- function(filename){
  #'@param filename Ensure that the directory from which you are loading the mogp has the 2 files "filename.RData" and "filename_mogp".
  #'@description Returns an ExUQmogp object: an R list with elements $mogp and $fitting.elements. The 2 elements must be loaded from 2 different files in this version of the code. These files will have been created with a call to save_ExUQmogp()
  ExUQmogpObj <- try(load(paste(filename,"RData",sep=".")),silent=TRUE)
  if(inherits(ExUQmogpObj,"try-error"))
    stop(paste("The file ", filename, ".RData does not exist. Ensure you specify filename so that both filename.RData and filename_mogp exist."))
  ExUQmogpObj <- eval(parse(text=ExUQmogpObj))
  if(!all(names(ExUQmogpObj) == c("mogp", "fitting.elements")))
    stop("Loaded object must be an Exeter UQ mogp object: a list with 2 elements: the mogp and the fitting elements.")
  ExUQmogpObj$mogp <- try(py_load_object(paste(filename,"mogp",sep="_")),silent=TRUE)
  if(inherits(ExUQmogpObj$mogp, "try-error"))
    stop(paste("The file ", filename, "_mogp does not exist. Ensure you specify filename so that both filename.RData and filename_mogp exist."))
  return(ExUQmogpObj)
}

virtual.LOO_MOGP <- function(mogp.emulator, lm.emulator) {
  #' @description Leave one out diagnostics computed using fast formula
  #' Function to compute the virtual Leave-One-Out formulas.
  #' @param mogp.emulator a single mogp emulator list
  #' @param lm.emulator linear model part of the emulator
  #' 
  #' @details Victoria to give the mathematical formulae being used. The lm.emulator is used to get the model matrix for the prior mean fit. This function may be added to mogp in the future. Need to check that mogp.emulator$nugget is always the right value even though a theta[last] value is calculated during modes such as "adaptive" and "fixed" that do not fit a nugget.
  
  #' @return  a random generation for the normal distribution with mean and standard deviation
  #' found using the virtual Leave-One-Out formulas.
  expectation <- c()
  variance <- c()
  
  H <- model.matrix(lm.emulator$linModel)
  n <- dim(H)[1]
  p <- dim(H)[2]
  Gamma2 <- mogp.emulator$get_K_matrix() + diag(mogp.emulator$nugget, n)
  Gamma2Inv <- solve(Gamma2)
  beta <- mogp.emulator$theta[1:p]
  MeanRes <- c(H%*%beta)
  
  Gamma2InvY <- Gamma2Inv %*%(mogp.emulator$targets - MeanRes)
  DiagGamma2Inv <- diag(Gamma2Inv)
  for(i in 1:n) {
    expectation[i] <- mogp.emulator$targets[i] - Gamma2InvY[i]/DiagGamma2Inv[i]
    variance[i] <- 1/DiagGamma2Inv[i]
  } 
  fit <- data.frame(cbind(expectation, expectation-2*sqrt(variance), 
                               expectation + 2*sqrt(variance)))
  names(fit) <- c('posterior mean', 'lower quantile', 'upper quantile')
  return(fit)
}

ValidationMOGP <- function(NewData, Emulators, which.emulator=1,
                           tData, ParamNames, 
                           Predictions=NULL, 
                           OriginalRanges = FALSE, 
                           RangeFile=NULL, Obs=NULL, 
                           ObsErr=NULL, ObsRange=FALSE) {
  #' @description Function to generate validation plots together with predictions for
  #' unseen data for a single emulator. Function also generates plots if predictions
  #' are provided.
  #' 
  #' @param NewData a data frame of unseen data that contains input parameters together
  #' with the model response of interest.
  #' @param Emulators a list of emulators fit by BuildNewEmulators
  #' @param which.emulator the index of the emulator you want predictions for.
  #' @param tData This is the data frame that we used to construct emulator. The format should be D+1+Q where Q is the number of targets each occupying one of the last Q columns, D is the number of inputs occupying the first D columns. The D+1th column should be a vector of random normal numbers (between -1 and 1) called "Noise". We use it in our LM code to stop our algorithms overfitting by selecting signals by chance. All variables should have names.
  #' @param Predictions A data frame with three columns with first column corresponding to posterior mean, 
  #' and second and third columns corresponding to the minus and plus two standard deviations.
  #'  If NULL then the predictions will be generated from emulator inside the function.
  #' @param ParamNames a vector of names for parameters.
  #' @param OriginalRanges. If TRUE, plots will be produced on the original parameter ranges
  #' Those ranges will be read from a file containing the parameter ranges and whether the 
  #' parameters are logged or not. Defaults to FALSE, where parameters will be plotted on [-1,1]
  #' @param RangeFile A .R file that will be sourced in order to determine the ranges of the
  #' parameters to be plotted and whether they are on a log scale or not. If NULL when OrignialRanges
  #' is called, a warning is thrown and the plot is given on [-1,1]
  #' @param Obs. The scalar value of the observations to be plotted as a dasked line if not NULL
  #' @param ObsErr. Observation error (scalar). If this is NULL when obs is not NULL, a warning is thrown.
  #' @param OriginalRanges. If TRUE, plots will be produced on the original parameter ranges
  #' Those ranges will be read from a file containing the parameter ranges and whether the 
  #' parameters are logged or not. Defaults to FALSE, where parameters will be plotted on [-1,1]
  #' @param RangeFile A .R file that will be sourced in order to determine the ranges of the
  #' parameters to be plotted and whether they are on a log scale or not. If NULL when OrignialRanges
  #' is called, a warning is thrown and the plot is given on [-1,1]
  #' @param ObsRange. Boolean to indicate if the plot window should have y ranges that include
  #' obs uncertainty (defaults to FALSE to facilitate emulator diagnostics)
  #' 
  #' @return a data frame with three columns, with first column corresponding to posterior mean, 
  #' and second and third columns corresponding to the minus and plus two standard deviations.
  #' The function also generates validation plots.
  # Description of a function goes here.
  lastCand <- which(names(tData)=="Noise")
  if(is.null(Predictions)) {
    predict.object <- Emulators$mogp$predict(as.matrix(NewData[, 1:(lastCand-1)]))
    if(Emulators$mogp$n==1) {
      fit.MOGP <- cbind(predict.object$mean[1, ], 
                        predict.object$mean[1, ]-2*sqrt(predict.object$unc[1, ]), 
                        predict.object$mean[1, ]+2*sqrt(predict.object$unc[1, ]))
      } else {
        fit.MOGP <- cbind(predict.object$mean[which.emulator, ], 
                          predict.object$mean[which.emulator, ]-2*sqrt(predict.object$unc[which.emulator, ]), 
                          predict.object$mean[which.emulator, ]+2*sqrt(predict.object$unc[which.emulator, ]))
      }
  } else {
    fit.MOGP <- Predictions
  }
  if(ObsRange){
    fit.MOGP <- rbind(fit.MOGP, c(Obs, Obs-2*ObsErr,Obs+2*ObsErr))
  }
  Design <- NewData[,Emulators$fitting.elements$ActiveIndices[[which.emulator]]]
  y <- NewData[, which(names(NewData) == Emulators$fitting.elements$lm.object[[which.emulator]]$ResponseString)]
  p <- length(ParamNames)
  if(p<2){
    par(mfrow = c(1, 1), mar=c(4, 4, 1, 1))
  }
  else if(p<3){
    par(mfrow = c(1, 2), mar=c(4, 4, 1, 1))
  }
  else if(p<4){
    par(mfrow = c(1, 3), mar=c(4, 4, 1, 1))
  }
  else if(p <5){
    par(mfrow = c(2, 2), mar=c(4, 4, 1, 1))
  }
  else if(p<7){
    par(mfrow = c(2, 3), mar=c(4, 4, 1, 1))
  }
  else if(p<10){
    par(mfrow = c(3, 3), mar=c(4, 4, 1, 1))
  }
  else if(p<13){
    par(mfrow = c(4, 3), mar=c(4, 4, 1, 1))
  }
  else if(p<=16){
    par(mfrow = c(4, 4), mar=c(4, 4, 1, 1))
  }
  if(OriginalRanges){
    if(is.null(RangeFile))
      stop("Cannot plot on original ranges as no RangeFile Specified")
    else{
      tRanFile <- try(source(RangeFile), silent=TRUE)
      if(inherits(tRanFile, "try-error"))
        stop("Invalid RangeFile given")
      if(!is.null(param.names)
         & !is.null(param.lows)
         & !is.null(param.highs)
         & !is.null(param.defaults)
         #& !is.null(which.logs)
      ){
        PlotOrder <- sapply(ParamNames, function(aName) which(param.names==aName))
        #TRY LAPPLY IF THE ABOVE FAILS NEEDING NUMERIC VECTORS NOT STRINGS
        #First cut design to just ParamNames in the order of ParamNames
        DesignOrder <- sapply(ParamNames, function(aName) which(colnames(Design)==aName))
        #Design order is a permutation with cut columns
        if(is.list(DesignOrder))
          DesignOrder <- unlist(DesignOrder)
        PermutedDesign <- Design[,DesignOrder]
        AllLogs <- rep(FALSE,length(param.names))
        AllLogs[which.logs] <- TRUE
        toInc <- which(names(PlotOrder)%in%colnames(PermutedDesign))
        param.names <- param.names[PlotOrder[toInc]]
        param.lows <- param.lows[PlotOrder[toInc]]
        param.highs <- param.highs[PlotOrder[toInc]]
        NewLogs <- AllLogs[PlotOrder[toInc]]
        which.logs <- which(NewLogs)
        Design <- DesignConvert(PermutedDesign, param.names = param.names, 
                                param.lows = param.lows, param.highs = param.highs, 
                                which.logs = which.logs)
      }
      else
        stop("Ranges file doesnt define the right variables")
    }
  }
  else{
    which.logs <- c()
  }
  tlogs <- rep("",p)
  tlogs[which.logs] <- "x"
  for(i in 1:p) {
    try(aplot <- ValidPlotNew(fit = fit.MOGP, x = Design[,i], 
                              y=y, ObsRange = ObsRange, 
                              main = "", cex.main=0.8,
                              xlab=ParamNames[i], log=tlogs[i]), silent=TRUE)
    if(!inherits(aplot, "try-error") & !is.null(Obs)){
      abline(h=Obs, lty=2, col=4)
      if(is.null(ObsErr))
        warning("The observations do not have 0 error. Please add ObsErr else the plot will be misleading")
      else{
        abline(h=Obs+2*ObsErr, col=4, lty=2)
        abline(h=Obs-2*ObsErr, col=4, lty=2)
      }  }
  }
  return(fit.MOGP)
}
LOO.plot <- function(Emulators, which.emulator=1, ParamNames,
                     OriginalRanges = FALSE, RangeFile=NULL, Obs=NULL, 
                     ObsErr=NULL, ObsRange=FALSE) {
  #Emulators, which.emulator=1, mogp.emulator, lm.emulator, Design, ParamNames,
   #                  OriginalRanges = FALSE, RangeFile=NULL, Obs=NULL, 
    #                 ObsErr=NULL, ObsRange=FALSE) {
  #' @description Function to generate Leave-One-Out plots for a whole design for a single emulator.
  #' 
  #' @param Emulators a list of emulators fit by BuildNewEmulators
  #' @param which.emulator the index of the emulator you want LOOs for.
  #' @param StanEmulator a GP emulator from EMULATE.gpstan function.
  #' @param ParamNames a vector of names for parameters.
  #' @param OriginalRanges. If TRUE, LOOs will be plotted on the original parameter ranges
  #' Those ranges will be read from a file containing the parameter ranges and whether the 
  #' parameters are logged or not. Defaults to FALSE, where parameters will be plotted on [-1,1]
  #' @param RangeFile A .R file that will be sourced in order to determine the ranges of the
  #' parameters to be plotted and whether they are on a log scale or not. If NULL when OrignialRanges
  #' is called, a warning is thrown and the plot is given on [-1,1]
  #' @param Obs. The scalar value of the observations to be plotted as a dasked line if not NULL
  #' @param ObsErr. Observation error (scalar). If this is NULL when obs is not NULL, a warning is thrown.
  #' @param OriginalRanges. If TRUE, LOOs will be plotted on the original parameter ranges
  #' Those ranges will be read from a file containing the parameter ranges and whether the 
  #' parameters are logged or not. Defaults to FALSE, where parameters will be plotted on [-1,1]
  #' @param RangeFile A .R file that will be sourced in order to determine the ranges of the
  #' parameters to be plotted and whether they are on a log scale or not. If NULL when OrignialRanges
  #' is called, a warning is thrown and the plot is given on [-1,1]
  #' @param ObsRange. Boolean to indicate if the plot window should have y ranges that include
  #' obs uncertainty (defaults to FALSE to facilitate emulator diagnostics)
  #' 
  #' @return a data frame with three columns, with first column corresponding to posterior mean, 
  #' and second and third columns corresponding to the minus and plus two standard deviations.
  #' The function also generates LOO validation plots.
  fit.loo <- virtual.LOO_MOGP(mogp.emulator=Emulators$mogp$emulators[[which.emulator]], lm.emulator=Emulators$fitting.elements$lm.object[[which.emulator]])
  if(ObsRange){
    fit.loo <- rbind(fit.loo, c(Obs, Obs-2*ObsErr,Obs+2*ObsErr))
  }
  Design <- Emulators$fitting.elements$Design[,Emulators$fitting.elements$ActiveIndices[[which.emulator]], drop = FALSE]
  p <- length(ParamNames)
  if(p<2){
    par(mfrow = c(1, 1), mar=c(4, 4, 1, 1))
  }
  else if(p<3){
    par(mfrow = c(1, 2), mar=c(4, 4, 1, 1))
  }
  else if(p<4){
    par(mfrow = c(1, 3), mar=c(4, 4, 1, 1))
  }
  else if(p <5){
    par(mfrow = c(2, 2), mar=c(4, 4, 1, 1))
  }
  else if(p<7){
    par(mfrow = c(2, 3), mar=c(4, 4, 1, 1))
  }
  else if(p<10){
    par(mfrow = c(3, 3), mar=c(4, 4, 1, 1))
  }
  else if(p<13){
    par(mfrow = c(4, 3), mar=c(4, 4, 1, 1))
  }
  else if(p<=16){
    par(mfrow = c(4, 4), mar=c(4, 4, 1, 1))
  }
  if(OriginalRanges){
    if(is.null(RangeFile))
      stop("Cannot plot on original ranges as no RangeFile Specified")
    else{
      tRanFile <- try(source(RangeFile), silent=TRUE)
      if(inherits(tRanFile, "try-error"))
        stop("Invalid RangeFile given")
      if(!is.null(param.names)
         & !is.null(param.lows)
         & !is.null(param.highs)
         & !is.null(param.defaults)
         #& !is.null(which.logs)
      ){
        PlotOrder <- sapply(ParamNames, function(aName) which(param.names==aName))
        #TRY LAPPLY IF THE ABOVE FAILS NEEDING NUMERIC VECTORS NOT STRINGS
        #First cut design to just ParamNames in the order of ParamNames
        DesignOrder <- sapply(ParamNames, function(aName) which(colnames(Design)==aName))
        #Design order is a permutation with cut columns
        if(is.list(DesignOrder))
          DesignOrder <- unlist(DesignOrder)
        PermutedDesign <- Design[,DesignOrder, drop=FALSE]
        AllLogs <- rep(FALSE,length(param.names))
        AllLogs[which.logs] <- TRUE
        toInc <- which(names(PlotOrder)%in%colnames(PermutedDesign))
        param.names <- param.names[PlotOrder[toInc]]
        param.lows <- param.lows[PlotOrder[toInc]]
        param.highs <- param.highs[PlotOrder[toInc]]
        NewLogs <- AllLogs[PlotOrder[toInc]]
        which.logs <- which(NewLogs)
        Design <- DesignConvert(PermutedDesign, param.names = param.names, 
                                   param.lows = param.lows, param.highs = param.highs, 
                                   which.logs = which.logs)
      }
      else
        stop("Ranges file doesnt define the right variables")
    }
  }
  else{
    which.logs <- c()
  }
  tlogs <- rep("",p)
  tlogs[which.logs] <- "x"
  for(i in 1:p) {
    try(aplot <- ValidPlotNew(fit = fit.loo, x = Design[,i], y=Emulators$mogp$emulators[[which.emulator]]$targets, 
                              ObsRange = ObsRange, main = "", cex.main=0.8,
                              xlab=ParamNames[i], log=tlogs[i]), silent=TRUE)
    if(!inherits(aplot, "try-error") & !is.null(Obs)){
      abline(h=Obs, lty=2, col=4)
      if(is.null(ObsErr))
        warning("The observations do not have 0 error. Please add ObsErr else the plot will be misleading")
      else{
        abline(h=Obs+2*ObsErr, col=4, lty=2)
        abline(h=Obs-2*ObsErr, col=4, lty=2)
      }      }
  }
  return(fit.loo)
}








ValidPlotNew <- function(fit, x, y, ObsRange=FALSE, ...){
  #' @param fit a data frame of emulator predictions. First column corresponds
  #' to the posterior mean, second and third columns correspond to the lower
  #' and upper quantiles (to be plotted) respectively.
  #' @param x values of the inputs to be plotted on the x axis
  #' @param y a vector of simulator evaluations at x
  #' @param ObsRange. Boolean to indicate if the last row of fit represents the
  #' observations and should not be plotted, but used to expand the ranges
  #' @param ... usual arguments to plotting functions. 
  inside <- c()
  outside <- c()
  for(i in 1:length(y)) {
    if(fit[i, 2] <= y[i] & y[i] <= fit[i, 3]) {
      inside <- c(inside,i)
    }
  }
  outside <- c(1:length(y))[-inside]
  tRange <- range(fit)
  if(ObsRange){
    fit <- fit[-nrow(fit),]
  }
  plot(x, fit[, 1], pch=20, ylim=tRange, ...)
  arrows(x, fit[, 2], x, fit[, 3],
         length=0.05, angle=90, code=3, col='black')
  points(x[inside], y[inside], pch=20, col='green')
  points(x[outside], y[outside], pch=20, col='red')
}

#Preprocessing for Optimal Basis emulation
ExtractCentredScaledDataAndBasis <- function(OriginalField, scaling=1){
  #MeanField is now the uncentered, unscaled field we want to emulate. The last step is to centre, scale and extract the basis
  EnsMean <- apply(OriginalField, MARGIN=1, mean)
  CentredField <- OriginalField
  for(i in 1:dim(OriginalField)[2]){
    CentredField[,i] <- CentredField[,i] - EnsMean
  }
  CentredField <- CentredField/scaling
  Basis <- svd(t(CentredField))$v
  return(list(tBasis=Basis, CentredField = CentredField, EnsembleMean = EnsMean, scaling=scaling))
}

##########################################################
#Design conversion functions (used in HighTune and for LOO plotting)
##########################################################
#All require specification of 
#param.lows (lower values of parameters)
#param.highs (max values of parameters)
#param.names (the names of the parameters)
#which.logs (pointers for those parameters that were designed on log scale)
#Mainly used in plotting

#Function to convert [-1,1] LHC to given scale
DesignConvert <- function(Xconts, param.names, param.lows, param.highs, which.logs=NULL){
  if(!(length(param.names)==length(param.lows)))
    stop("specify as many parameter names as parameter ranges")
  else if(!(length(param.highs)==length(param.lows)))
    stop("Mismatch in min and max values")
  conversion <- function(anX,lows,highs){
    ((anX+1)/2)*(highs-lows) +lows
  }
  param.lows.log <- param.lows
  param.highs.log <- param.highs
  param.lows.log[which.logs] <- log10(param.lows[which.logs])
  param.highs.log[which.logs] <- log10(param.highs[which.logs])
  tX <- sapply(1:length(param.lows), function(i) conversion(Xconts[,i],param.lows.log[i],param.highs.log[i]))
  tX[,which.logs] <- 10^tX[,which.logs]
  tX <- as.data.frame(tX)
  names(tX) <- param.names
  tX
}

#3. Function to convert [above scale to [-1, 1]
#Probably need to pass parameter ranges to be safe
DesignantiConvert <- function (Xconts){
  anticonversion <- function(newX,lows,highs){
    2*((newX-lows)/(highs-lows))-1
  }
  param.lows.log <- param.lows
  param.highs.log <- param.highs
  param.lows.log[which.logs] <- log10(param.lows[which.logs])
  param.highs.log[which.logs] <- log10(param.highs[which.logs])
  Xconts[,which.logs] <- log10(Xconts[,which.logs])
  tX <- sapply(1:length(param.lows), function(i) anticonversion(Xconts[,i],param.lows.log[i],param.highs.log[i]))
  tX <- as.data.frame(tX)
  names(tX) <- param.names
  tX
}
DesignantiConvert1D <- function (Xconts){
  anticonversion <- function(newX,lows,highs){
    2*((newX-lows)/(highs-lows))-1
  }
  param.lows.log <- param.lows
  param.highs.log <- param.highs
  param.lows.log[which.logs] <- log10(param.lows[which.logs])
  param.highs.log[which.logs] <- log10(param.highs[which.logs])
  Xconts[which.logs] <- log10(Xconts[which.logs])
  tX <- sapply(1:length(param.lows), function(i) anticonversion(Xconts[i],param.lows.log[i],param.highs.log[i]))
  tX <- as.data.frame(tX)
  tX
}


