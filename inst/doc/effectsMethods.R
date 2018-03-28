## ----include=FALSE-------------------------------------------------------
library(knitr)
opts_chunk$set(
tidy=FALSE,fig.width=5,fig.height=5,cache=FALSE
)

## ----echo=FALSE, results='hide', include=FALSE---------------------------
#options(continue="+    ", prompt="R> ", width=76)
options(show.signif.stars=FALSE)
options(scipen=3)

## ----include=FALSE-------------------------------------------------------
library(knitr)
opts_chunk$set(
tidy=FALSE,fig.width=5,fig.height=5,cache=FALSE,comment=NA, prompt=TRUE
)
render_sweave()

## ----echo=FALSE, results='hide', include=FALSE----------------------------
options(continue="+    ", prompt="R> ", width=76)
options(show.signif.stars=FALSE)
options(scipen=3)

## ----eval=FALSE-----------------------------------------------------------
#  Effect.gls <- function(focal.predictors, mod, ...){
#    args <- list(
#      type = "glm",
#      call = mod$call,
#      formula = formula(mod),
#      family = family(mod),
#      link = NULL,
#      coefficients = coef(mod),
#      vcov = as.matrix(vcov(mod)))
#    Effect.default(focal.predictors, mod, ..., sources=args)
#  }

## ----eval=FALSE-----------------------------------------------------------
#  Effect.merMod <- function(focal.predictors, mod, ..., KR=FALSE){
#    if (KR && !requireNamespace("pbkrtest", quietly=TRUE)){
#      KR <- FALSE
#      warning("pbkrtest is not available, KR set to FALSE")}
#    fam <- family(mod)
#    args <- list(
#      call = mod@call,
#      coefficients = lme4::fixef(mod),
#      vcov = if (fam == "gaussian" && fam$link == "identity" && KR)
#        as.matrix(pbkrtest::vcovAdj(mod)) else as.matrix(vcov(mod)))
#    Effect.default(focal.predictors, mod, ..., sources=args)
#  }

## ----eval=FALSE-----------------------------------------------------------
#  Effect.rlmerMod <- function(focal.predictors, mod, ...){
#    args <- list(
#      coefficients = lme4::fixef(mod))
#    Effect.default(focal.predictors, mod, ..., sources=args)
#  }

## ----eval=FALSE-----------------------------------------------------------
#  Effect.betareg <- function(focal.predictors, mod, ...){
#    coef <- mod$coefficients$mean
#    vco <- vcov(mod)[1:length(coef), 1:length(coef)]
#    args <- list(
#      call = mod$call,
#      formula = formula(gy_logit),
#      coefficients = coef,
#      link = mod$link$mean,
#      vcov = vco)
#    Effect.default(focal.predictors, mod, ..., sources=args)
#  }

## -------------------------------------------------------------------------
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

## ----eval=FALSE-----------------------------------------------------------
#  Effect.clm2 <- function(focal.predictors, mod, ...){
#    if (requireNamespace("MASS", quietly=TRUE)){
#      polr <- MASS::polr}
#    if(is.null(mod$Hessian)){
#      message("\nRe-fitting to get Hessian\n")
#      mod <- update(mod, Hess=TRUE)}
#    if(mod$link != "logistic") stop("Effects only supports the logit link")
#    if(mod$threshold != "flexible") stop("Effects only supports the flexible threshold")
#    numTheta <- length(mod$Theta)
#    numBeta <- length(mod$beta)
#    or <- c( (numTheta+1):(numTheta + numBeta), 1:(numTheta))
#    args <- list(
#      type = "polr",
#      formula = mod$call$location,
#      coefficients = mod$beta,
#      vcov = as.matrix(vcov(mod)[or, or]))
#    Effect.default(focal.predictors, mod, ..., sources=args)
#  }

## ----eval=FALSE-----------------------------------------------------------
#  Effect.clmm <- function(focal.predictors, mod, ...){
#    if (requireNamespace("MASS", quietly=TRUE)){
#      polr <- MASS::polr}
#    if(is.null(mod$Hessian)){
#      message("\nRe-fitting to get Hessian\n")
#      mod <- update(mod, Hess=TRUE)}
#    if(mod$link != "logit") stop("Only the logistic link is supported by Effects")
#    if(mod$threshold != "flexible") stop("Only threshold='flexible supported by Effects")
#    numTheta <- length(mod$Theta)
#    numBeta <- length(mod$beta)
#    or <- c( (numTheta+1):(numTheta + numBeta), 1:(numTheta))
#    skip <- length(unique(model.frame(mod)[,1])) - 1
#    vcov <- matrix(NA, nrow=numBeta + skip, ncol=numBeta + skip)
#    sel <- rownames(vcov(mod)) %in% names(mod$beta)
#    vcov[1:numBeta, 1:numBeta] <- vcov(mod)[sel, sel]
#    args <- list(
#      type = "polr",
#      formula = fixFormula(as.formula(mod$formula)),
#      coefficients = mod$beta,
#      vcov = as.matrix(vcov))
#    Effect.default(focal.predictors, mod, ..., sources=args)
#  }

