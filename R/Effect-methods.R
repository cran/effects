# 12/11/2017:  This file contains all the Effect methods that call Effect.default.  S. Weisberg
#              Excluded are Effect.lm, Effect.polr, and Effect.multinom, and for now Effect.svyglm.

# new lme method
Effect.lme <- function(focal.predictors, mod, ...){
  args <- list(
    call = mod$call,
    formula = mod$call$fixed,
    coefficients = mod$coefficients$fixed,
    vcov = as.matrix(vcov(mod)))
  Effect.default(focal.predictors, mod, ..., sources=args)
}

# new gls method
Effect.gls <- function(focal.predictors, mod, ...){
  args <- list(
    call = mod$call,
    formula = formula(mod),
    coefficients = coef(mod),
    vcov = as.matrix(vcov(mod)))
  Effect.default(focal.predictors, mod, ..., sources=args)
}

# new merMod
Effect.merMod <- function(focal.predictors, mod, ..., KR=FALSE){
  if (KR && !requireNamespace("pbkrtest", quietly=TRUE)){
    KR <- FALSE
    warning("pbkrtest is not available, KR set to FALSE")}
  fam <- family(mod)
  args <- list(
    call = mod@call,
    coefficients = lme4::fixef(mod),
    vcov = if (fam == "gaussian" && fam$link == "identity" && KR)
      as.matrix(pbkrtest::vcovAdj(mod)) else as.matrix(vcov(mod)))
  Effect.default(focal.predictors, mod, ..., sources=args)
}

# rlmer in robustlmm package
Effect.rlmerMod <- function(focal.predictors, mod, ...){
  args <- list(
    coefficients = lme4::fixef(mod))
  Effect.default(focal.predictors, mod, ..., sources=args)
}

# clm in the ordinal package
Effect.clm <- function(focal.predictors, mod, ...){
  if (requireNamespace("MASS", quietly=TRUE)){
    polr <- MASS::polr}
  if(mod$link != "logit") stop("Effects only supports the logit link")
  if(mod$threshold != "flexible") stop("Effects only supports the flexible threshold")
  if(is.null(mod$Hessian)){
    message("\nRe-fitting to get Hessian\n")
    mod <- update(mod, Hess=TRUE)}
  numTheta <- length(mod$Theta)
  numBeta <- length(mod$beta)
  or <- c( (numTheta+1):(numTheta + numBeta), 1:(numTheta))
  args <- list(
    type = "polr",
    coefficients = mod$beta,
    vcov = as.matrix(vcov(mod)[or, or]))
  Effect.default(focal.predictors, mod, ..., sources=args)
}

# clm2
Effect.clm2 <- function(focal.predictors, mod, ...){
  if (requireNamespace("MASS", quietly=TRUE)){
      polr <- MASS::polr}
  if(is.null(mod$Hessian)){
     message("\nRe-fitting to get Hessian\n")
     mod <- update(mod, Hess=TRUE)}
  if(mod$link != "logistic") stop("Effects only supports the logit link")
  if(mod$threshold != "flexible") stop("Effects only supports the flexible threshold")
  numTheta <- length(mod$Theta)
  numBeta <- length(mod$beta)
  or <- c( (numTheta+1):(numTheta + numBeta), 1:(numTheta))
  args <- list(
    type = "polr",
    formula = mod$call$location,
    coefficients = mod$beta,
    vcov = as.matrix(vcov(mod)[or, or]))
  Effect.default(focal.predictors, mod, ..., sources=args)
}

#clmm in ordinal package
Effect.clmm <- function(focal.predictors, mod, ...){
  if (requireNamespace("MASS", quietly=TRUE)){
    polr <- MASS::polr}
  if(is.null(mod$Hessian)){
    message("\nRe-fitting to get Hessian\n")
    mod <- update(mod, Hess=TRUE)}
  if(mod$link != "logit") stop("Only the logistic link is supported by Effects")
  if(mod$threshold != "flexible") stop("Only threshold='flexible supported by Effects")
  numTheta <- length(mod$Theta)
  numBeta <- length(mod$beta)
  or <- c( (numTheta+1):(numTheta + numBeta), 1:(numTheta))
  skip <- length(unique(model.frame(mod)[,1])) - 1
  vcov <- matrix(NA, nrow=numBeta + skip, ncol=numBeta + skip)
  sel <- rownames(vcov(mod)) %in% names(mod$beta)
  vcov[1:numBeta, 1:numBeta] <- vcov(mod)[sel, sel]
  args <- list(
    type = "polr",
    formula = fixFormula(as.formula(mod$formula)),
    coefficients = mod$beta,
    vcov = as.matrix(vcov))
  Effect.default(focal.predictors, mod, ..., sources=args)
}

# betareg from the betareg package
Effect.betareg <- function(focal.predictors, mod, ...){
  coef <- mod$coefficients$mean
  vco <- vcov(mod)[1:length(coef), 1:length(coef)]
  args <- list(
    call = mod$call,
    formula = formula(mod),
    coefficients = coef,
    link = mod$link$mean,
    vcov = vco)
  Effect.default(focal.predictors, mod, ..., sources=args)
}


