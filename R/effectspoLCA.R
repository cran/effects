# 2013-07-31:  extend effects to poLCA objects.  S. Weisberg
# 2013-10-15: removed effect.poLCA. J. Fox
# 2018-11-19: added focal.levels argument to predictorEffect() and predictorEffects() methods. J. Fox

#The next two functions should be exported to the namespace

allEffects.poLCA <- function(mod, ...){
  allEffects(poLCA.to.fake(mod), ...)
}

Effect.poLCA <- function(focal.predictors, mod, ...) {
  result <- Effect(focal.predictors, poLCA.to.fake(mod), ..., sources=list(type="multinom"))
  result$formula <- as.formula(formula(mod))
  result
}

predictorEffects.poLCA <- function(mod, predictors = ~.,focal.levels=50, xlevels=5,  ...){
  predictorEffects(poLCA.to.fake(mod), predictors=predictors, focal.levels=focal.levels, xlevels=xlevels, ...)
}

predictorEffect.poLCA <- function(predictor, mod, focal.levels=50, xlevels=5, ...){
  predictorEffect(predictor, poLCA.to.fake(mod), focal.levels=focal.levels, xlevels=xlevels, ...)
}

# this function makes a 'fake' multinom object or 'glm' object so
# effect.multinom  or effect.glm can be used.
# effect.multinom requires at least 3 classes, so if classes=2 use
# effect.glm
poLCA.to.fake <- function(mod) { 
  if (requireNamespace("nnet", quietly=TRUE)){
    multinom <- nnet::multinom}
  dta <- eval(mod$call$data)
  form <- as.formula(eval(mod$call$formula))
  # find the missing data:
  omit <- attr(model.frame(form, dta), "na.action")
  if(length(omit) == 0) dta$.class <- factor(mod$predclass) else{
    dta$.class <- rep(NA, dim(dta)[1])
    dta$.class[-omit] <- mod$predclass
    dta$.class <- factor(dta$.class)
  }
  # end of missing data correction
  formula1 <- update(form, .class ~ .)
  if(length(mod$P) == 2L){
    mod1 <- glm(formula1, family=binomial, data=dta)
    mod1$call$data <- dta
    mod1$call$formula <- formula1
    mod1$coef <- mod$coeff[, 1]
    mod1$vcov <- mod$coeff.V
    class(mod1) <- c("fakeglm", class(mod1)) }
  else {
    mod1 <- multinom(formula1, dta, Hess=TRUE, trace=FALSE, maxit=1)
    mod1$call$data <- dta
    mod1$call$formula <- formula1
    mod1$coeff <- mod$coeff
    mod1$coeff.V <- mod$coeff.V
    class(mod1) <- c("fakemultinom", class(mod1))
  }
  coef.fakemultinom <- function(mod){
    coef <- t(mod$coeff)
    dimnames(coef) <- list(mod$lab[-1L], mod$vcoefnames)
    coef
  }
  vcov.fakemultinom <- function(mod){mod$coeff.V}
  mod1
}



