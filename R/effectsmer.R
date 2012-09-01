# effect.mer and effect.lme built from effect.lm by S. Weisberg 29 June 2011
# last modified 2012-03-08 to require() lme4 or nlme. J. Fox


# the function lm.wfit fit gets the hessian wrong for mer's.  Get the variance
# from the vcov method applied to the mer object.

# 'fixmod' is a copy of the 'nobars' function in the lme4 package, 
# renamed so it doesn't cause any conflicts.  This is a utility function
# that should not be exported

fixmod <- function (term) 
{
    if (!("|" %in% all.names(term))) 
        return(term)
    if (is.call(term) && term[[1]] == as.name("|")) 
        return(NULL)
    if (length(term) == 2) {
        nb <- fixmod(term[[2]])
        if (is.null(nb)) 
            return(NULL)
        term[[2]] <- nb
        return(term)
    }
    nb2 <- fixmod(term[[2]])
    nb3 <- fixmod(term[[3]])
    if (is.null(nb2)) 
        return(nb3)
    if (is.null(nb3)) 
        return(nb2)
    term[[2]] <- nb2
    term[[3]] <- nb3
    term
}

# mer.to.glm evaluates a 'glm' model that is as similar to a given 'mer'
# model as follows.  It is of class c("fakeglm", "glm", "lm")
# several items are added to the created objects. Do not export
mer.to.glm <- function(mod) {
    cl <- mod@call
    if(cl[[1]] =="nlmer") stop("effects package does not support 'nlmer' objects")
    m <- match(c("formula", "family", "data", "weights", "subset", 
                 "na.action", "start", "offset",  
                 "model", "contrasts"), names(cl), 0L)
    cl <- cl[c(1L, m)]
    cl[[1L]] <- as.name("glm")
    cl$formula <- fixmod(as.formula(cl$formula))
    mod2 <- eval(cl)
    mod2$coefficients <- mod@fixef
    mod2$vcov <- as.matrix(vcov(mod))
    mod2$linear.predictors <- model.matrix(mod2) %*% mod2$coefficients
    mod2$fitted.values <- mod2$family$linkinv(mod2$linear.predictors)
    mod2$weights <- as.vector(with(mod2,
          prior.weights * (family$mu.eta(linear.predictors)^2 /
                           family$variance(fitted.values))))
    mod2$residuals <- with(mod2,
          prior.weights * (y - fitted.values)/weights )
    class(mod2) <- c("fakeglm", class(mod2))
    mod2
    }
# lme.to.glm evaluates a 'glm' model that is as similar to a given 'lme'
# model, in the same pattern as mer.to.glm.  This could be speeded up
# slightly by using 'lm' rather than 'glm' but I use 'glm' to parallel
# mer.to.glm more closely.  The differences are:  (1) match fewer args
# in the call; (2) different def of mod2$coefficients; no other 
# changes
lme.to.glm <- function(mod) {
    cl <- mod$call
    cl$formula <- cl$fixed
    m <- match(c("formula", "data", "weights", "subset", 
                 "na.action",  "contrasts"), names(cl), 0L)
    cl <- cl[c(1L, m)]
    cl[[1L]] <- as.name("glm")
    mod2 <- eval(cl)
    mod2$coefficients <- mod$coefficients$fixed
    mod2$vcov <- as.matrix(vcov(mod))
    mod2$linear.predictors <- model.matrix(mod2) %*% mod2$coefficients
    mod2$fitted.values <- mod2$family$linkinv(mod2$linear.predictors)
    mod2$weights <- as.vector(with(mod2,
          prior.weights * (family$mu.eta(linear.predictors)^2 /
                           family$variance(fitted.values))))
    mod2$residuals <- with(mod2,
          prior.weights * (y - fitted.values)/weights )
    class(mod2) <- c("fakeglm", class(mod2))
    mod2
    }

#method for 'fakeglm' objects. Do not export   
vcov.fakeglm <- function(object, ...) object$vcov

#The next four functions should be exported
effect.mer <- function(term, mod, ...) {
	if (!require(lme4)) stop("the lme4 package is not installed")
    result <- effect(term, mer.to.glm(mod), ...)
    result$formula <- as.formula(formula(mod))
    result
    }
    
allEffects.mer <- function(mod, ...){
  	if (!require(lme4)) stop("the lme4 package is not installed")
  	allEffects(mer.to.glm(mod), ...)
}
 
allEffects.lme <- function(mod, ...){
	if (!require(nlme)) stop("the nlme package is not installed")
  	allEffects(lme.to.glm(mod), ...)
}
  
effect.lme <- function(term, mod, ...) {
	if (!require(nlme)) stop("the nlme package is not installed")
    result <- effect(term, lme.to.glm(mod), ...)
    result$formula <- as.formula(formula(mod))
    result
    }
   
