# Effect generic and methods
# John Fox and Sanford Weisberg
# 2012-12-21: Allow for empty cells in factor interactions, S. Weisberg
# 2012-03-05: Added .merMod method for development version of lme4, J. Fox
# 2012-04-06: Added support for lme4.0, J. Fox
# 2013-07-15:  Changed default xlevels and default.levels
# 2013-10-15: Added Effect.default(). J. Fox
# 2013-10-22: fixed bug in Effect.lm() when na.action=na.exclude. J. Fox
# 2013-10-29: code to handle "valid" NAs in factors. J. Fox
# 2013-11-06: fixed bug in Effect.multinom() in construction of effect object
# 2014-03-13: modified Effect.lm() to compute partial residuals. J. Fox
# 2014-05-06: fixed bug in Effect.gls() when cor or var structure depends on variables in the data set. J. Fox
# 2014-08-02: added vcov.=vcov argument to allow other methods of estimating var(coef.estimates)
# 2014-09-25: added KR argument to Effect.mer() and Effect.merMod(). J. Fox
# 2014-12-07: don't assume that pbkrtest is installed. J. Fox
# 2015-03-25: added "family" element to eff objects returned by Effect.lm(). J. Fox
# 2016-02-16: fixed problem in handling terms like polynomials for non-focal predictors. J. Fox
# 2016-03-01: recoded calculation of partial residuals. J. Fox
# 2016-07-19: added checkFormula(). J. Fox
# 2017-08-18: removed default.levels argument. J. Fox
# 2017-08-26: introduced confint list argument, including Scheffe intervals. J. Fox
# 2017-08-29: reintroduce legacy se and confidence.level arguments.
# 2017-09-07: added Effect.svyglm()
# 2017-09-14: no partial residuals for Effect.svyglm()
# 2017-11-03: correct handling of rank deficient models, now using `estimability` package
# 2017-11-22: modified checkFormula to work with clm2 models that don't have a 'formula' argument
# 2017-12-10: Effect.default. Effect.mer, .merMod, .lme, gls have been replaced to use the default.
# 2018-01-22: allow given.values="equal" or given.values="default" 
# 2018-01-25: substitute se for confint arg; make confint a legacy arg
# 2018-05-06: allow for complete=FALSE arg in potential calls to vcov.lm() and vcov.glm.
# 2018-05-13: allow partial residuals to be computed when the x.var is a factor.
# 2018-06-05: Effect.default now makes sure family$aic is 
#             set, for use with non-standard families.
# 2018-06-05: A test has been added to Effect.default to chech if family$variance
#             has one parameter.  If not, the function is stopped and an error is
#             returned.
# 2018-06-12: Fixed bug with vcov in Effect.default
# 2018-06-20: Added a check to Effect.default to handle family args that
#             are character or an unevaluated function
# 2018-10-01: Avoid warnings when testing given.values == "equal" or "default".
# 2018-10-08: transformation argument changed to legacy
# 2018-10-08: new returned value 'link' = family(mod) 
# 2019-04-20: made Effect.default() more robust in fitting fake glm by setting epsilon=Inf.
# 2019-04-20: fixed bug in .set.given.equal() in tests for model class.
# 2019-07-05: clm, clm2 and clmm were not passing threshholds to the fake polr object, now corrected.
# 2019-09-04: handle xlevels=n argument correctly

### Non-exported function added 2018-01-22 to generalize given.values to allow for "equal" weighting of factor levels for non-focal predictors.
.set.given.equal <- function(m){
  if(inherits(m, "lm") & !("(Intercept)" %in% names(coef(m))))
      stop("Seting given.vales='equal' requires an intercept in the model formula")
  terms <- terms(m)
  classes <- attr(terms, "dataClasses")
  response <- attr(terms, "response")
  classes <- classes[-response]
  factors <- names(classes)[classes=="factor"]
  out <- NULL
  for (f in factors){
    form <- as.formula(paste( "~", f, collapse=""))
    .m0 <- if(inherits(m, "glm")) 
              {update(m, form, control=glm.control(epsilon=Inf, maxit=1))} else {
    if(inherits(m, "polr"))
              {update(m, form, control=list(maxit=1))} else {
    if(inherits(m, "multinom"))    
              {update(m, form, maxit=0, trace=FALSE)} else
        update(m, form)}}
    names <- colnames(model.matrix(.m0))[-1]
    vals <- rep(1/(length(names)+1), length(names))
    names(vals) <- names
    out <- c(out, vals)
  }
  out
}
### end of non-exported function

checkFormula <- function(object){
# clm2 does not have a formula,
  if(inherits(object, "clm2")) formula <- function(x) x$call$location
  if (!inherits(object, "formula")){
    object <- formula(object)
  }
  formula <- as.character(object)
  rhs <- formula[length(formula)]
  res <- regexpr("as.factor\\(|factor\\(|as.ordered\\(|ordered\\(|as.numeric\\(|as.integer\\(",
                 rhs)
  res == -1 || attr(res, "match.length") == 0
}

Effect <- function(focal.predictors, mod, ...){
  if (!checkFormula(mod)) stop("model formula should not contain calls to",
                               "\n  factor(), as.factor(), ordered(), as.ordered(),",
                               " as.numeric(), or as.integer();",
                               "\n  see 'Warnings and Limitations' in ?Effect")
  UseMethod("Effect", mod)
}

# 2017-12-04 new Effect.default that actually works
# 2017-12-07 added Effects.lme, .mer, gls that work

Effect.default <- function(focal.predictors, mod, ..., sources=NULL){ 
# get formula from sources if present else from mod
  formula <- fixFormula(
        if(is.null(sources$formula)) formula(mod) else sources$formula)
# the next line returns the formula if focal.predictors is null
  if(is.null(focal.predictors)) return(formula)
# get the call
  cl <- if(is.null(sources$call)) {if(isS4(mod)) 
        mod@call else mod$call} else sources$call
# insert formula into the call
  cl$formula <- formula
# set type == 'glm' unless it is set in sources
  type <- if(is.null(sources$type)) "glm" else sources$type
# glm family from sources if set, else set fam to NULL
  if(!is.null(sources$family)){
    fam <- sources$family
    fam$aic <- function(...) NULL
    # check to be sure the variance function in the family has one argument only,
    # otherwise this method won't work
    if(!is.null(fam$variance)){
      if(length(formals(fam$variance)) > 1)
        stop("Effect plots are not implemented for families with more than
             one parameter in the variance function (e.g., negitave binomials).")}
  } else {fam <- NULL}
# get the coefficient estimates and vcov from sources if present
  coefficients <- if(is.null(sources$coefficients)) 
    coef(mod) else sources$coefficients
# added 7/5/2019, next line, for models that use polr (e.g, clm, clm2)
  zeta <- if(is.null(sources$zeta)) NULL else sources$zeta
  vcov <- if(is.null(sources$vcov)) 
    as.matrix(vcov(mod, complete=TRUE)) else sources$vcov
# end reading sources
# set control parameters: suggested by Nate TeGrotenhuis
  cl$control <- switch(type,
          glm = glm.control(epsilon=Inf, maxit=1),
          polr = list(maxit=1),
          multinom = c(maxit=1))
  cl$method <- sources$method # NULL except for type=="polr"
  .m <- switch(type,
               glm=match(c("formula", "data", "contrasts",  "subset",
                "control", "offset"), names(cl), 0L),
               polr=match(c("formula", "data", "contrasts",  "subset",
                "control", "method"), names(cl), 0L),
               multinom=match(c("formula", "data", "contrasts",  "subset",
                "family", "maxit", "offset"), names(cl), 0L))
  cl <- cl[c(1L, .m)]
  if(!is.null(fam)) cl$family <- fam
  if (is.character(cl$family)) 
    cl$family <- get(cl$family, mode = "function", envir = parent.frame())
  if (is.function(cl$family)) 
    cl$family <- family()
  cl[[1L]] <- as.name(type)
# The following eval creates on object of class glm, polr or multinom.  
# These are crated to avoid writing an Effects method for every type of model.  
# The only information used from this "fake" object are the coefficients and 
# the variance-covariance matrix, and these are copied from the original 
# object so Effects plots the right things.
  mod2 <- eval(cl)
  mod2$coefficients <- coefficients
  mod2$vcov <- vcov
  if(!is.null(zeta)) mod2$zeta <- zeta # added 7/5/2019
  if(type == "glm"){
       mod2$weights <- as.vector(with(mod2,
                        prior.weights * (family$mu.eta(linear.predictors)^2 /
                                        family$variance(fitted.values))))}
  class(mod2) <- c("fakeeffmod", class(mod2))
  Effect(focal.predictors, mod2, ...)  # call the glm/polr/multinom method
}

vcov.fakeeffmod <- function(object, ...) object$vcov

## This function removes terms with "|" or "||" in the formula, assumking these
## correspond to random effects.
fixFormula <- function (term)
{
  if (!("|" %in% all.names(term)) && !("||" %in% all.names(term)))
    return(term)
  if ((is.call(term) && term[[1]] == as.name("|")) ||
      (is.call(term) && term[[1]] == as.name("||")))
    return(NULL)
  if (length(term) == 2) {
    nb <- fixFormula(term[[2]])
    if (is.null(nb))
      return(NULL)
    term[[2]] <- nb
    return(term)
  }
  nb2 <- fixFormula(term[[2]])
  nb3 <- fixFormula(term[[3]])
  if (is.null(nb2))
    return(nb3)
  if (is.null(nb3))
    return(nb2)
  term[[2]] <- nb2
  term[[3]] <- nb3
  term
}

Effect.lm <- function(focal.predictors, mod, xlevels=list(), fixed.predictors,
        vcov. = vcov, se=TRUE,
        residuals=FALSE, quantiles=seq(0.2, 0.8, by=0.2),
        x.var=NULL,  ...,
        #legacy arguments:
        given.values, typical, offset, confint, confidence.level, 
        partial.residuals, transformation){ 
  
  if (is.numeric(xlevels)){
    if (length(xlevels) > 1 || round(xlevels != xlevels)) stop("xlevels must be a single whole number or a list")
    form <- Effect.default(NULL, mod) #returns the fixed-effects formula
    terms <- attr(terms(form), "term.labels")
    predictors <- all.vars(parse(text=terms))
    xlevs <- list()
    for (pred in predictors){
      xlevs[[pred]] <- xlevels
    }
    xlevels <- xlevs
  }
  
  if (!missing(partial.residuals)) residuals <- partial.residuals
  partial.residuals <- residuals
  if (missing(transformation)) 
    transformation <- list(link = family(mod)$linkfun, 
                           inverse = family(mod)$linkinv)
  if (missing(fixed.predictors)) fixed.predictors <- NULL
  fixed.predictors <- applyDefaults(fixed.predictors,
                                    list(given.values=NULL, typical=mean,
                                         apply.typical.to.factors=FALSE, offset=mean),
                                    arg="fixed.predictors")
  if (missing(given.values)) given.values <- fixed.predictors$given.values
# new 1/22/18 to allow for automatical equal weighting of factor levels
  if(!is.null(given.values)){
   if (given.values[1] == "default") given.values <- NULL
   if (given.values[1] == "equal") given.values <- .set.given.equal(mod)}
# end new code
  if (missing(typical)) typical <- fixed.predictors$typical
  if (missing(offset)) offset <- fixed.predictors$offset
  apply.typical.to.factors <- fixed.predictors$apply.typical.to.factors
  if (!missing(confint)) se <- confint
  confint <- applyDefaults(se, list(compute=TRUE, level=.95, type="pointwise"),
                           onFALSE=list(compute=FALSE, level=.95, type="pointwise"),
                           arg="se")
  se <- confint$compute
  if (missing(confidence.level)) confidence.level <- confint$level
  confidence.type <- match.arg(confint$type, c("pointwise", "Scheffe", "scheffe"))
  default.levels <- NULL # just for backwards compatibility
  data <- if (partial.residuals){
    all.vars <- all.vars(formula(mod))
    expand.model.frame(mod, all.vars)[, all.vars]
  }
  else NULL
  if (!is.null(given.values) && !all(which <- names(given.values) %in% names(coef(mod))))
    stop("given.values (", names(given.values[!which]), ") not in the model")
  off <- if (is.numeric(offset) && length(offset) == 1) offset
  else if (is.function(offset)) {
    mod.off <- model.offset(model.frame(mod))
    if (is.null(mod.off)) 0 else offset(mod.off)
  }
  else stop("offset must be a function or a number")
  formula.rhs <- formula(mod)[[3]]
  if (!missing(x.var)){
    if (!is.numeric(x.var)) {
      x.var.name <- x.var
      x.var <- which(x.var == focal.predictors)
    }
    if (length(x.var) == 0) stop("'", x.var.name, "' is not among the focal predictors")
    if (length(x.var) > 1) stop("x.var argument must be of length 1")
  }
  model.components <- Analyze.model(focal.predictors, mod, xlevels, default.levels, formula.rhs,
                                    partial.residuals=partial.residuals, quantiles=quantiles, x.var=x.var, data=data, typical=typical)
  excluded.predictors <- model.components$excluded.predictors
  predict.data <- model.components$predict.data
  predict.data.all.rounded <- predict.data.all <- if (partial.residuals) na.omit(data[, all.vars(formula(mod))]) else NULL
  factor.levels <- model.components$factor.levels
  factor.cols <- model.components$factor.cols
  n.focal <- model.components$n.focal
  x <- model.components$x
  X.mod <- model.components$X.mod
  cnames <- model.components$cnames
  X <- model.components$X
  x.var <- model.components$x.var
  formula.rhs <- formula(mod)[c(1, 3)]
  Terms <- delete.response(terms(mod))
  mf <- model.frame(Terms, predict.data, xlev = factor.levels, na.action=NULL)
  mod.matrix <- model.matrix(formula.rhs, data = mf, contrasts.arg = mod$contrasts)
  if (is.null(x.var)) partial.residuals <- FALSE
  factors <- sapply(predict.data, is.factor)
  if (partial.residuals){
    for (predictor in focal.predictors[-x.var]){
      if (!factors[predictor]){
        values <- unique(predict.data[, predictor])
        predict.data.all.rounded[, predictor] <- values[apply(outer(predict.data.all[, predictor], values, function(x, y) (x - y)^2), 1, which.min)]
      }
    }
  }
  mod.matrix.all <- model.matrix(mod)
  wts <- weights(mod)
  if (is.null(wts))
    wts <- rep(1, length(residuals(mod)))
  mod.matrix <- Fixup.model.matrix(mod, mod.matrix, mod.matrix.all,
                                   X.mod, factor.cols, cnames, focal.predictors, 
                                   excluded.predictors, typical, given.values, 
                                   apply.typical.to.factors) 
# 11/3/2017.  Check to see if the model is full rank
  # Compute a basis for the null space, using estimibility package
  null.basis <- estimability::nonest.basis(mod)  # returns basis for null space
  # check to see if each row of mod.matrix is estimable
  is.estimable <- estimability::is.estble(mod.matrix, null.basis) # TRUE if effect is estimable else FALSE
  # substitute 0 for NA in coef vector and compute effects
  scoef <- ifelse(is.na(mod$coefficients), 0L, mod$coefficients)
  effect <- off + mod.matrix %*% scoef
  effect[!is.estimable] <- NA  # set all non-estimable effects to NA
# end estimability check
  if (partial.residuals){
    res <- na.omit(residuals(mod, type="working"))
    fitted <- na.omit(if (inherits(mod, "glm")) predict(mod, type="link") else predict(mod))
    partial.residuals.range <- range(fitted + res)
  }
  else {
    res <- partial.residuals.range <- NULL
  }
  result <- list(term = paste(focal.predictors, collapse="*"),
                 formula = formula(mod), response = response.name(mod),
                 variables = x, fit = effect, x = predict.data[, 1:n.focal, drop=FALSE],
                 x.all=predict.data.all.rounded[, focal.predictors, drop=FALSE],
                 model.matrix = mod.matrix,
                 data = X,
                 discrepancy = 0, offset=off,
                 residuals=res, partial.residuals.range=partial.residuals.range,
                 x.var=x.var)
  if (se) {
    if (any(family(mod)$family == c("binomial", "poisson"))) {
      z <- if (confidence.type == "pointwise") {
        qnorm(1 - (1 - confidence.level)/2)
      } else {
        p <- length(na.omit(coef(mod)))
        scheffe(confidence.level, p)
      }
    }
    else {
      z <- if (confidence.type == "pointwise") {
        qt(1 - (1 - confidence.level)/2, df = mod$df.residual)
      } else {
        p <- length(na.omit(coef(mod)))
        scheffe(confidence.level, p, mod$df.residual)
      }
    }
    V <- vcov.(mod, complete=FALSE)
    mmat <- mod.matrix[, !is.na(mod$coefficients)] # remove non-cols with NA coeffs
    eff.vcov <- mmat %*% V %*% t(mmat)
    rownames(eff.vcov) <- colnames(eff.vcov) <- NULL
    var <- diag(eff.vcov)
    result$vcov <- eff.vcov
    result$se <- sqrt(var)
    result$se[!is.estimable] <- NA
    result$lower <- effect - z * result$se
    result$upper <- effect + z * result$se
    result$confidence.level <- confidence.level
  }
  if (is.null(transformation$link) && is.null(transformation$inverse)) {
    transformation$link <- I
    transformation$inverse <- I
  }
  result$transformation <- transformation
  result$family <- family(mod)$family
# 2018-10-08 result$family kept to work with legacy code  
  result$link <- family(mod) 
  class(result) <- "eff"
  result
}

Effect.multinom <- function(focal.predictors, mod,
                            xlevels=list(), fixed.predictors,
                            vcov. = vcov, se=TRUE, ...,
                            #legacy arguments:
                            confint, confidence.level, given.values, typical){
  
  if (is.numeric(xlevels)){
    if (length(xlevels) > 1 || round(xlevels != xlevels)) stop("xlevels must be a single whole number or a list")
    form <- Effect.default(NULL, mod) #returns the fixed-effects formula
    terms <- attr(terms(form), "term.labels")
    predictors <- all.vars(parse(text=terms))
    xlevs <- list()
    for (pred in predictors){
      xlevs[[pred]] <- xlevels
    }
    xlevels <- xlevs
  }
  
  if (missing(fixed.predictors)) fixed.predictors <- NULL
  fixed.predictors <- applyDefaults(fixed.predictors,
                                    list(given.values=NULL, typical=mean),
                                    arg="fixed.predictors")
  if (missing(given.values)) given.values <- fixed.predictors$given.values
  # new 1/22/18 to allow for automatical equal weighting of factor levels
  if(!is.null(given.values)){
    if (given.values[1] == "default") given.values <- NULL
    if (given.values[1] == "equal") given.values <- .set.given.equal(mod)}
  # end new code
  # end new code
  if (missing(typical)) typical <- fixed.predictors$typical
  if (!missing(confint)) se <- confint
  confint <- applyDefaults(se, list(compute=TRUE, level=.95, type="pointwise"),
                           onFALSE=list(compute=FALSE, level=.95, type="pointwise"),
                           arg="se")
  se <- confint$compute
  if (missing(confidence.level)) confidence.level <- confint$level
  confidence.type <- match.arg(confint$type, c("pointwise", "Scheffe", "scheffe"))
  default.levels <- NULL # just for backwards compatibility
  if (length(mod$lev) < 3) stop("effects for multinomial logit model only available for response levels > 2")
  if (missing(given.values)) given.values <- NULL
  else if (!all(which <- colnames(given.values) %in% names(coef(mod))))
    stop("given.values (", colnames(given.values[!which]),") not in the model")
  formula.rhs <- formula(mod)[c(1, 3)]
  model.components <- Analyze.model(focal.predictors, mod, xlevels, default.levels, formula.rhs, typical=typical)
  excluded.predictors <- model.components$excluded.predictors
  predict.data <- model.components$predict.data
  factor.levels <- model.components$factor.levels
  factor.cols <- model.components$factor.cols
  #    n.focal <- model.components$n.focal
  x <- model.components$x
  X.mod <- model.components$X.mod
  cnames <- model.components$cnames
  X <- model.components$X
  formula.rhs <- formula(mod)[c(1, 3)]
  Terms <- delete.response(terms(mod))
  mf <- model.frame(Terms, predict.data, xlev = factor.levels)
  mod.matrix <- model.matrix(formula.rhs, data = mf, contrasts.arg = mod$contrasts)
  X0 <- Fixup.model.matrix(mod, mod.matrix, model.matrix(mod),
                           X.mod, factor.cols, cnames, focal.predictors, excluded.predictors, typical, given.values)
  resp.names <- make.names(mod$lev, unique=TRUE)
  resp.names <- c(resp.names[-1], resp.names[1]) # make the last level the reference level
  B <- t(coef(mod))
  V <- vcov.(mod)
  m <- ncol(B) + 1
  p <- nrow(B)
  r <- p*(m - 1)
  n <- nrow(X0)
  P <- Logit <- matrix(0, n, m)
  colnames(P) <-  paste("prob.", resp.names, sep="")
  colnames(Logit) <-  paste("logit.", resp.names, sep="")
  if (se){
    z <- if (confidence.type == "pointwise") {
      qnorm(1 - (1 - confidence.level)/2)
    } else {
      scheffe(confidence.level, p)
    }
    Lower.P <- Upper.P <- Lower.logit <- Upper.logit <- SE.P <- SE.logit <- matrix(0, n, m)
    colnames(Lower.logit) <-  paste("L.logit.", resp.names, sep="")
    colnames(Upper.logit) <-  paste("U.logit.", resp.names, sep="")
    colnames(Lower.P) <-  paste("L.prob.", resp.names, sep="")
    colnames(Upper.P) <-  paste("U.prob.", resp.names, sep="")
    colnames(SE.P) <-  paste("se.prob.", resp.names, sep="")
    colnames(SE.logit) <-  paste("se.logit.", resp.names, sep="")
  }
  for (i in 1:n){
    res <- eff.mul(X0[i,], B, se, m, p, r, V) # compute effects
    #        P[i,] <- prob <- res$p # fitted probabilities
    P[i,] <- res$p # fitted probabilities
    Logit[i,] <- logit <- res$logits # fitted logits
    if (se){
      #            SE.P[i,] <- se.p <- res$std.err.p # std. errors of fitted probs
      SE.P[i,] <- res$std.err.p # std. errors of fitted probs
      SE.logit[i,] <- se.logit <- res$std.error.logits # std. errors of logits
      Lower.P[i,] <- logit2p(logit - z*se.logit)
      Upper.P[i,] <- logit2p(logit + z*se.logit)
      Lower.logit[i,] <- logit - z*se.logit
      Upper.logit[i,] <- logit + z*se.logit
    }
  }
  resp.levs <- c(m, 1:(m-1)) # restore the order of the levels
  P <- P[, resp.levs]
  Logit <- Logit[, resp.levs]
  if (se){
    Lower.P <- Lower.P[, resp.levs]
    Upper.P <- Upper.P[, resp.levs]
    Lower.logit <- Lower.logit[, resp.levs]
    Upper.logit <- Upper.logit[, resp.levs]
    SE.P <- SE.P[, resp.levs]
    SE.logit <- SE.logit[, resp.levs]
  }
  result <- list(term=paste(focal.predictors, collapse="*"), formula=formula(mod), response=response.name(mod),
                 y.levels=mod$lev, variables=x, x=predict.data[, focal.predictors, drop=FALSE],
                 model.matrix=X0, data=X, discrepancy=0, model="multinom",
                 prob=P, logit=Logit)
  if (se) result <- c(result, list(se.prob=SE.P, se.logit=SE.logit,
                                   lower.logit=Lower.logit, upper.logit=Upper.logit,
                                   lower.prob=Lower.P, upper.prob=Upper.P,
                                   confidence.level=confidence.level))
  # find empty cells, if any, and correct
## 11/3/17:  The code until the next comment is surely incorrect, but
## generally harmless.  One must learn if the notion of estimablilty applied
## to multinomial models and figure out the right thing to do
  whichFact <- unlist(lapply(result$variables, function(x) x$is.factor))
  zeroes <- NULL
  if(sum(whichFact) > 1){
    nameFact <- names(whichFact)[whichFact]
    counts <- xtabs(as.formula( paste("~", paste(nameFact, collapse="+"))),
                    model.frame(mod))
    zeroes <- which(counts == 0)
  }
  if(length(zeroes) > 0){
    levs <- expand.grid(lapply(result$variables, function(x) x$levels))
    good <- rep(TRUE, dim(levs)[1])
    for(z in zeroes){
      good <- good &
        apply(levs, 1, function(x) !all(x == levs[z, whichFact]))
    }
    result$prob[!good, ] <- NA
    result$logit[!good, ] <- NA
    if (se){
      result$se.prob[!good, ] <- NA
      result$se.logit[!good, ] <- NA
      result$lower.prob[!good, ] <- NA
      result$upper.prob[!good, ] <- NA
    }
  }
## End of unnecessary code
  class(result) <-'effpoly'
  result
}

Effect.polr <- function(focal.predictors, mod,
                        xlevels=list(), fixed.predictors,
                        vcov.=vcov, se=TRUE, latent=FALSE, ...,
                        #legacy arguments:
                        confint, confidence.level, given.values, typical){
  
  if (is.numeric(xlevels)){
    if (length(xlevels) > 1 || round(xlevels != xlevels)) stop("xlevels must be a single whole number or a list")
    form <- Effect.default(NULL, mod) #returns the fixed-effects formula
    terms <- attr(terms(form), "term.labels")
    predictors <- all.vars(parse(text=terms))
    xlevs <- list()
    for (pred in predictors){
      xlevs[[pred]] <- xlevels
    }
    xlevels <- xlevs
  }
  
  if (missing(fixed.predictors)) fixed.predictors <- NULL
  fixed.predictors <- applyDefaults(fixed.predictors,
                                    list(given.values=NULL, typical=mean),
                                    arg="fixed.predictors")
  if (missing(given.values)) given.values <- fixed.predictors$given.values
  # new 1/22/18 to allow for automatical equal weighting of factor levels
  # new 1/22/18 to allow for automatical equal weighting of factor levels
  if(!is.null(given.values)){
    if (given.values[1] == "default") given.values <- NULL
    if (given.values[1] == "equal") given.values <- .set.given.equal(mod)}
  # end new code
  if (missing(typical)) typical <- fixed.predictors$typical
  if (!missing(confint)) se <- confint
  confint <- applyDefaults(se, list(compute=TRUE, level=.95, type="pointwise"),
                           onFALSE=list(compute=FALSE, level=.95, type="pointwise"),
                           arg="se")
  se <- confint$compute
  if (missing(confidence.level)) confidence.level <- confint$level
  confidence.type <- match.arg(confint$type, c("pointwise", "Scheffe", "scheffe"))
  default.levels <- NULL # just for backwards compatibility
  if (mod$method != "logistic") stop('method argument to polr must be "logistic"')
  if (missing(given.values)) given.values <- NULL
  else if (!all(which <- names(given.values) %in% names(coef(mod))))
    stop("given.values (", names(given.values[!which]),") not in the model")
  formula.rhs <- formula(mod)[c(1, 3)]
  model.components <- Analyze.model(focal.predictors, mod, xlevels, default.levels, formula.rhs, typical=typical)
  excluded.predictors <- model.components$excluded.predictors
  predict.data <- model.components$predict.data
  factor.levels <- model.components$factor.levels
  factor.cols <- model.components$factor.cols
  #    n.focal <- model.components$n.focal
  x <- model.components$x
  X.mod <- model.components$X.mod
  cnames <- model.components$cnames
  X <- model.components$X
  Terms <- delete.response(terms(mod))
  mf <- model.frame(Terms, predict.data, xlev = factor.levels, na.action=NULL)
  mod.matrix <- model.matrix(formula.rhs, data = mf, contrasts.arg = mod$contrasts)
  X0 <- Fixup.model.matrix(mod, mod.matrix, model.matrix(mod),
                           X.mod, factor.cols, cnames, focal.predictors, excluded.predictors, typical, given.values)
  resp.names <- make.names(mod$lev, unique=TRUE)
  X0 <- X0[,-1, drop=FALSE]
  b <- coef(mod)
  p <- length(b)  # corresponds to p - 1 in the text
  alpha <- - mod$zeta  # intercepts are negatives of thresholds
  z <- if (confidence.type == "pointwise") {
    qnorm(1 - (1 - confidence.level)/2)
  } else {
    scheffe(confidence.level, p + length(alpha))
  }
  result <- list(term=paste(focal.predictors, collapse="*"), formula=formula(mod), response=response.name(mod),
                 y.levels=mod$lev, variables=x,
                 x=predict.data[, focal.predictors, drop=FALSE],
                 model.matrix=X0, data=X, discrepancy=0, model="polr")
  if (latent){
    res <- eff.latent(X0, b, vcov.(mod)[1:p, 1:p], se)
    result$fit <- res$fit
    if (se){
      result$se <- res$se
      result$lower <- result$fit - z*result$se
      result$upper <- result$fit + z*result$se
      result$confidence.level <- confidence.level
    }
    transformation <- list()
    transformation$link <- I
    transformation$inverse <- I
    result$transformation <- transformation
    result$thresholds <- -alpha
    class(result) <- c("efflatent", "eff")
    return(result)
  }
  m <- length(alpha) + 1
  r <- m + p - 1
  indices <- c((p+1):r, 1:p)
  V <- vcov.(mod)[indices, indices]
  for (j in 1:(m-1)){  # fix up the signs of the covariances
    V[j,] <- -V[j,]  #  for the intercepts
    V[,j] <- -V[,j]}
  n <- nrow(X0)
  P <- Logit <- matrix(0, n, m)
  colnames(P) <-  paste("prob.", resp.names, sep="")
  colnames(Logit) <-  paste("logit.", resp.names, sep="")
  if (se){
    Lower.logit <- Upper.logit <- Lower.P <- Upper.P <- SE.P <- SE.Logit <- matrix(0, n, m)
    colnames(Lower.logit) <-  paste("L.logit.", resp.names, sep="")
    colnames(Upper.logit) <-  paste("U.logit.", resp.names, sep="")
    colnames(Lower.P) <-  paste("L.prob.", resp.names, sep="")
    colnames(Upper.P) <-  paste("U.prob.", resp.names, sep="")
    colnames(SE.P) <-  paste("se.prob.", resp.names, sep="")
    colnames(SE.Logit) <-  paste("se.logit.", resp.names, sep="")
  }
  for (i in 1:n){
    res <- eff.polr(X0[i,], b, alpha, V, m, r, se) # compute effects
    P[i,] <- res$p # fitted probabilities
    Logit[i,] <- logit <- res$logits # fitted logits
    if (se){
      SE.P[i,] <- res$std.err.p # std. errors of fitted probs
      SE.Logit[i,] <- se.logit <- res$std.error.logits # std. errors of logits
      Lower.P[i,] <- logit2p(logit - z*se.logit)
      Upper.P[i,] <- logit2p(logit + z*se.logit)
      Lower.logit[i,] <- logit - z*se.logit
      Upper.logit[i,] <- logit + z*se.logit
    }
  }
  result$prob <- P
  result$logit <- Logit
  if (se) result <- c(result,
                      list(se.prob=SE.P, se.logit=SE.Logit,
                           lower.logit=Lower.logit, upper.logit=Upper.logit,
                           lower.prob=Lower.P, upper.prob=Upper.P,
                           confidence.level=confidence.level))
  class(result) <-'effpoly'
  result
}

# svyglm
Effect.svyglm <- function(focal.predictors, mod, fixed.predictors, ...){
  Svymean <- function(x){
    svymean(x, design=mod$survey.design)
  }
  ellipses.list <- list(...)
  if ((!is.null(ellipses.list$residuals) && !isFALSE(residuals)) || 
      (!is.null(ellipses.list$partial.residuals) && !isFALSE(ellipses.list$partial.residuals))){
    stop("partial residuals are not available for svyglm models")
  }
  if (missing(fixed.predictors)) fixed.predictors <- NULL
  fixed.predictors <- applyDefaults(fixed.predictors,
                                    list(given.values=NULL, typical=Svymean,
                                         apply.typical.to.factors=TRUE, offset=Svymean),
                                    arg="fixed.predictors")
  typical <- fixed.predictors$typical
  apply.typical.to.factors <- fixed.predictors$apply.typical.to.factors
  offset <- fixed.predictors$offset
  mod$call <- list(mod$call, data=mod$data)
  Effect.lm(focal.predictors, mod, typical=typical,
            apply.typical.to.factors=apply.typical.to.factors, offset=offset, ...)
}

