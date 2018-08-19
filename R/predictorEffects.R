# last update:  2017-07-16
# 2017-08-14 fixed bug in plot.predictoreff on passing 'multiline' to lines list
# 2017-08-30 for compatibility with other effect plots, default
#            is now multiline=FALSE
# 2017-11-09 fixed bug in setting the class for multinom models, and possibly others
# 2017-11-17 added methods for clm, clm2, clmm in the file effectsclmm.R
# 2017-12-08 modified predictorEffect.default and predictorEffects.default for compatibility to Effect.default
# 2018-01-09 fixed bug in predictorEffects.default with log() in a formula.
# 2018-01-24 fixed bug with minus sign in a formula predictorEffects.default
# 2018-05-14 predictorEffect.default() calls Effect() with x.var=1
# 2018-06-07 predictorEffects now works with offsets.
# 2018-08-09 removed explicit 'xlevels' argument from predictorEffects, so the argument is correctly passed with ...


# removed xlevels argument 8/9/18
predictorEffect <- function(predictor, mod, ...){  
  UseMethod("predictorEffect", mod)
}

# removed xlevels argument 8/9/18
predictorEffect.svyglm <- function(predictor, mod, ...){
  mod$call <- list(mod$call, data=mod$data)
  NextMethod(object=mod, ...)
}

#simplified 12/10/17
# removed xlevels argument 8/9/18
predictorEffect.default <- function(predictor, mod, ...){
  form <- Effect.default(NULL, mod) #returns the fixed-effects formula
  all.vars <- all.vars(parse(text=form))
  # find the right effect to use
  terms <- attr(terms(form), "term.labels")
  # get the predictor names:
  predictors <- all.vars(parse(text=terms))
  sel <- which(predictors == predictor)
  if(length(sel) != 1) stop("First argument must be the quoted name of one predictor in the formula")
  # create correspondence table
  decode <- function(name) all.vars(parse(text=unlist(strsplit(name, ":"))))
  tab <- rep(FALSE, length(terms))
  for(j in 1:length(terms)){if(predictor %in% decode(terms[j])) tab[j] <- TRUE}
  ans <- unlist(strsplit(paste(terms[tab], collapse=":"), ":"))
  ans <- unique(all.vars(parse(text=ans)))
  ans <- unique(c(predictor, ans)) # guarantees focal predictor is first
  args <- names(list(...))
  result <- if ("x.var" %in% args) Effect(ans, mod, ...) else Effect(ans, mod, x.var=1, ...)
  class(result) <- c("predictoreff", class(result))
  result
}

predictorEffects <- function(mod, predictors, ...){
  UseMethod("predictorEffects", mod)
}


# rewritten, simplified, 12/08/17, bug in formulas fixed 01/24/2018
predictorEffects.default <- function(mod, predictors = ~ ., ...) {
# The next function removes offset(s) from a formula, used for mform and cform
  no.offset <- function(x, preserve = NULL) {
    k <- 0
    proc <- function(x) {
      if (length(x) == 1) return(x)
      if (x[[1]] == as.name("offset") && !((k<<-k+1) %in% preserve)) return(x[[1]])
      replace(x, -1, lapply(x[-1], proc))
    }
    update(proc(x), . ~ . - offset)}
  mform <- no.offset(Effect.default(NULL, mod))  # returns the fixed-effect formula for any method
  cform <- if(is.character(predictors)) 
    as.formula(paste("~", paste(predictors, collapse="+"))) else
      predictors
  cform <- update(as.formula(paste(". ~", 
                paste(all.vars(formula(mform)[[3]]), collapse="+"))), 
                cform)
  cform <- no.offset(cform)
  mvars <- all.vars(mform[[3]])
  cvars <- all.vars(cform[[3]])
# check that 'cvars' is a subset of 'mvars'. If so apply predictorEffect
  if(!all(cvars %in% mvars))
    stop("argument 'predictors' not a subset of the predictors in the formula") else {
    result <- list()
    for(p in cvars) result[[p]] <- predictorEffect(p, mod, ...)
    }
  class(result) <- c("predictorefflist", "list")
  result
}

# plot methods

plot.predictoreff <- function(x, x.var,
                              main = paste(names(x$variables)[1], "predictor effect plot"), ...){
  if(missing(x.var)) x.var <- names(x$variables)[1]
  NextMethod(x, x.var=x.var, main=main, ...)
}


plot.predictorefflist <- function(x, selection, rows, cols, ask=FALSE, graphics=TRUE,
                                  lattice, ...){
  # Next line added 8/23/17 along with lattice, also lattice arg above
  lattice <- if(missing(lattice)) list() else lattice
  if (!missing(selection)){
    if (is.character(selection)) selection <- gsub(" ", "", selection)
    return(plot(x[[selection]], ...))
  }
  effects <- gsub(":", "*", names(x))
  if (ask){
    repeat {
      selection <- menu(effects, graphics=graphics, title="Select Term to Plot")
      if (selection == 0) break
      else print(plot(x[[selection]], ...))
    }
  }
  else {
    neffects <- length(x)
    mfrow <- mfrow(neffects)
    if (missing(rows) || missing(cols)){
      rows <- mfrow[1]
      cols <- mfrow[2]
    }
    for (i in 1:rows) {
      for (j in 1:cols){
        if ((i-1)*cols + j > neffects) break
        more <- !((i-1)*cols + j == neffects)
        lattice[["array"]] <- list(row=i, col=j, nrow=rows, ncol=cols, more=more)
        print(plot(x[[(i-1)*cols + j]], lattice=lattice,
                   ...))
      }
    }
  }
}


# print and summary methods

print.predictorefflist <- function(x, ...){
  for (eff in x){
    print(eff, ...)
  }
  invisible(x)
}

print.predictoreff <- function(x, ...){
  cat("\n", names(x$variables)[1], "predictor effect\n")
  NextMethod()
}

summary.predictorefflist <- function(object, ...){
  for (eff in object){
    cat("\n", names(eff$variables)[1], "predictor effect\n")
    print(summary(eff, ...))
  }
}
