# last update:  2017-07-16
# 2017-08-14 fixed bug in plot.predictoreff on passing 'multiline' to lines list
# 2017-08-30 for compatibility with other effect plots, default 
#            is now multiline=FALSE

predictorEffect <- function(predictor, mod, xlevels, ...){
  UseMethod("predictorEffect", mod)
}

predictorEffect.svyglm <- function(predictor, mod, xlevels, ...){
  mod$call <- list(mod$call, data=mod$data)
  NextMethod(object=mod)
}

predictorEffect.default <- function(predictor, mod, xlevels=list(), ...){
  all.vars <- all.vars(formula(mod))
  data <- na.omit(expand.model.frame(mod, all.vars)[, all.vars])
  if (is.null(xlevels[[predictor]]) && is.numeric(data[[predictor]])){
    xlevels[[predictor]] <- quantile(data[[predictor]], seq(0.01, 0.99, by=0.02))
  } 
  # find the right effect to use
  terms <- attr(terms(mod), "term.labels")
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
  result <- Effect(ans, mod, xlevels=xlevels, ...)
  class(result) <- c("predictoreff", "eff")
  result
}

predictorEffects <- function(mod, predictors = ~ ., ...){
  # convert `predictors` arg to a list of predictors
  vterms <- if(is.character(predictors)) paste("~",predictors) else predictors
  vform <- update(formula(mod), vterms)
  vlabels <- attr(terms(vform), "term.labels")
  vpred <- all.vars(parse(text=vlabels))
  # get list of predictors from the model
  mlabels <- attr(attr(model.frame(mod), "terms"), "term.labels")
  mpred <- all.vars(parse(text=mlabels))
  # check that 'vpred' is a subset of 'mpred'. If so apply predictorEffect
  if(!all(vpred %in% mpred)) stop("argument 'predictors' not a subset of the predictors") else {
    result <- list()
    for(p in vpred) result[[p]] <- predictorEffect(p, mod, ...)
  }
  class(result) <- 'predictorefflist'
  result
}

# plot methods

plot.predictoreff <- function(x, x.var, 
                              main = paste(names(x$variables)[1], "predictor effect plot"), ...){
  if(missing(x.var)) x.var <- names(x$variables)[1]
  NextMethod(x, x.var=x.var, main=main, ...)
}

# This next function differs for plot.efflist only by changing the title on each plot
plot.predictorefflist <- function(x, selection, rows, cols, ask=FALSE, graphics=TRUE, 
                                  lattice, ...){
  lattice <- if(missing(lattice)) list() else lattice
  if(length(x) == 1) plot(x[[1]],  ...) else {
    if (!missing(selection)){
      if (is.character(selection)) selection <- gsub(" ", "", selection)
      return(plot(x[[selection]], ...))
    }
    effects <- gsub(":", "*", names(x))
    if (ask){
      repeat {
        selection <- menu(effects, graphics=graphics, title="Select Term to Plot")
        if (selection == 0) break
        else print(plot(x[[selection]], z.var=names(x)[selection], ...))
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
                     x.var=names(x)[(i-1)*cols + j],
                     main=paste(names(x)[(i-1)*cols + j], "predictor effect plot"), ...))
        }
      }
    }
  }}

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
