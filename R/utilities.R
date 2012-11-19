# utilities and common functions for effects package
# John Fox, Jangman Hong, and Sanford Weisberg
#  last modified 2012-11-19 by J. Fox

if (getRversion() >= "2.15.1") globalVariables("wt")


has.intercept <- function(model, ...) any(names(coefficients(model))=="(Intercept)")

term.names <- function (model, ...) {
	term.names <- gsub(" ", "", labels(terms(model)))
	if (has.intercept(model)) c("(Intercept)", term.names)
	else term.names
}

response.name <- function (model, ...) deparse(attr(terms(model), "variables")[[2]])

mfrow <- function(n, max.plots=0){
	# number of rows and columns for array of n plots
	if (max.plots != 0 & n > max.plots)
		stop(paste("number of plots =",n," exceeds maximum =", max.plots))
	rows <- round(sqrt(n))
	cols <- ceiling(n/rows)
	c(rows, cols)
}

expand.model.frame <- function (model, extras, envir = environment(formula(model)),
	na.expand = FALSE){  # modified version of R base function
	f <- formula(model)
	data <- eval(model$call$data, envir)
	ff <- foo ~ bar + baz
	if (is.call(extras)) 
		gg <- extras
	else gg <- parse(text = paste("~", paste(extras, collapse = "+")))[[1]]
	ff[[2]] <- f[[2]]
	ff[[3]][[2]] <- f[[3]]
	ff[[3]][[3]] <- gg[[2]]
	if (!na.expand) {
		naa <- model$call$na.action
		subset <- model$call$subset
		rval <- if (is.null(data)) eval(call("model.frame", ff, # modified
						subset = subset, na.action = naa), envir)           #  lines
			else eval(call("model.frame", ff, data = data,          #
						subset = subset, na.action = naa), envir)           #
	}
	else {
		subset <- model$call$subset
		rval <- eval(call("model.frame", ff, data = data, subset = subset, 
				na.action = I), envir)
		oldmf <- model.frame(model)
		keep <- match(rownames(oldmf), rownames(rval))
		rval <- rval[keep, ]
		class(rval) <- "data.frame"
	}
	return(rval)
}

is.relative <- function(term1, term2, factors) {
	all(!(factors[,term1]&(!factors[,term2])))
}

ancestors <- function(term, mod,...){
	names <- term.names(mod)
	if (has.intercept(mod)) names <- names[-1]
	if(length(names) == 1) return(NULL)
	which.term <- which(term == names)
	if (length(which.term) == 0){
		factors <- attr(terms(...), "factors")
		rownames(factors) <- gsub(" ", "", rownames(factors))
		colnames(factors) <- gsub(" ", "", colnames(factors))
		result<-(1:length(names))[sapply(names,
				function(term2) is.relative(term2, term, factors))]
		if (0 ==  length(result)) which.term else result
	}
	else {
		factors <- attr(terms(mod), "factors")     
		rownames(factors) <- gsub(" ", "", rownames(factors))
		colnames(factors) <- gsub(" ", "", colnames(factors))   
		result<-(1:length(names))[-which.term][sapply(names[-which.term],
				function(term2) is.relative(term2, term, factors))]
		if (0 ==  length(result)) which.term else result
	}
}

descendants <- function(term, mod, ...){
	names <- term.names(mod)
	if (has.intercept(mod)) names <- names[-1]
	if(length(names)==1) return(NULL)
	which.term <- which(term == names)
	if (length(which.term) == 0){
		factors <- attr(terms(...), "factors")
		rownames(factors) <- gsub(" ", "", rownames(factors))
		colnames(factors) <- gsub(" ", "", colnames(factors))
		(1:length(names))[sapply(names,
				function(term2) is.relative(term, term2, factors))]
	}
	else {
		factors <- attr(terms(mod), "factors")
		rownames(factors) <- gsub(" ", "", rownames(factors))
		colnames(factors) <- gsub(" ", "", colnames(factors))
		(1:length(names))[-which.term][sapply(names[-which.term],
				function(term2) is.relative(term, term2, factors))]
	}
}

first.order.ancestors <- function(term, mod,...){
	ancestors <- ancestors(term, mod, ...)
	ancestors[attr(terms(mod), 'order')[ancestors] == 1]
}

is.high.order.term <- function(term, mod,...){
	0 == length(descendants(term, mod, ...))
}

subscripts <- function(index, dims){
	subs <- function(dims, index){
		dim <- length(dims)
		if (dim == 0) return(NULL)
		cum <- c(1,cumprod(dims))[dim]
		i <- index %/% cum
		if (index %% cum != 0) i <- i + 1
		c(i, subs(dims[-dim], index - (i - 1)*cum))
	}
	rev(subs(dims, index))
}

#matrix.to.df <- function(matrix){
#	on.exit(options(warn = opt[[1]]))
#	opt <- options(warn = -1)
#	ncol <- ncol(matrix)
#	colnames <- colnames(matrix)
#	result <- list()
#	for (j in 1:ncol){
#		numbers <- as.numeric(matrix[,j])
#		result[[colnames[j]]] <-
#			if(all(is.na(numbers))) matrix[,j] else numbers
#	}
#	as.data.frame(result)
#}

matrix.to.df <- function(matrix, colclasses){
	on.exit(options(warn = opt[[1]]))
	opt <- options(warn = -1)
	ncol <- ncol(matrix)
	colnames <- colnames(matrix)
	colclasses[sapply(colclasses, function(x) "integer" %in% x)] <- "numeric"
	result <- vector(mode="list", length=ncol)
	names(result) <- colnames
	for (j in 1:ncol){
		result[[j]] <- matrix[, j]
		class <- colclasses[[colnames[j]]]
		result[[colnames[j]]] <- if ("numeric" %in% class) {
					decChar <- getOption('OutDec')
					if (decChar == '.') as.numeric(result[[colnames[j]]])
					else as.numeric(gsub(decChar, '.', matrix[,j]))
				}
				else if ("ordered" %in% class) ordered(result[[colnames[j]]])
				else if ("factor" %in% class) factor(result[[colnames[j]]]) 
				else result[[colnames[j]]]
	}
	as.data.frame(result)
}

strangers <- function(term, mod,...){
	names <- term.names(mod)
	if (has.intercept(mod)) names <- names[-1]
	self <- which(names == term)
	ancestors <- ancestors(term, mod, ...)
	descendants <- descendants(term, mod, ...)
	sort(setdiff(1:ncol(attr(terms(mod), "factors")),
			union(union(ancestors, descendants), self)))
}

analyze.model <- function(term, mod, xlevels, default.levels){
    if ((!is.null(mod$na.action)) && class(mod$na.action) == "exclude") 
        class(mod$na.action) <- "omit"
    term <- gsub(" ", "", gsub("\\*", ":", term))
    intercept <- has.intercept(mod)
    terms <- term.names(mod)
    if (intercept) terms <- terms[-1]
    which.term <- which(term == terms)
    mod.aug<- list()
    if (length(which.term) == 0){
        warning(paste(term,"does not appear in the model"))
        mod.aug <- update(formula(mod), eval(parse(text=paste(". ~ . +", term))))
    }
    if (!is.high.order.term(term, mod, mod.aug))
        warning(paste(term, 'is not a high-order term in the model'))
    basic.vars <- first.order.ancestors(term, mod, mod.aug)
    all.vars <- (1:nrow(attr(terms(mod), 'factors')))[
        0 != apply(attr(terms(mod), 'factors'), 1, sum) ]
    if (intercept) all.vars <- all.vars - 1
    if (inherits(mod, "multinom")) all.vars <- all.vars - 1
    if (inherits(mod, "polr")) all.vars <- all.vars - 1
    
    excluded.vars <- setdiff(all.vars, basic.vars) 
    if (length(terms) == 1) {
        all.vars <- basic.vars <- tail(all.vars(formula(mod)), 1)
        excluded.vars <- numeric()
    }
    else {
        all.vars <- all.vars(as.formula(paste ("~", paste(terms[all.vars], collapse="+"))))
        basic.vars <- all.vars(as.formula(paste ("~", paste(terms[basic.vars], collapse="+"))))
    }
    excluded.vars <- if (length(excluded.vars) > 0) 
        all.vars(as.formula(paste ("~", paste(terms[excluded.vars], collapse="+"))))
    else NULL
    X.mod <- model.matrix(mod)
    cnames <- colnames(X.mod)
    factor.cols <- rep(FALSE, length(cnames))
    names(factor.cols) <- cnames
    X <- model.frame(mod)
    for (name in all.vars){
        if (is.factor(X[[name]])) factor.cols[grep(paste("^", name, sep=""), cnames)] <- TRUE
    }
    factor.cols[grep(":", cnames)] <- FALSE   
    X <- na.omit(expand.model.frame(mod, all.vars))
    x<-list()
    factor.levels <- list()
    for (name in basic.vars){
        levels <- mod$xlevels[[name]]
        fac <- !is.null(levels)
        if (!fac) {
            levels <- if (is.null(xlevels[[name]]))
                seq(min(X[, name]), max(X[,name]), length=default.levels)
            else xlevels[[name]]
        }
        else factor.levels[[name]] <- levels
        x[[name]] <- list(name=name, is.factor=fac, levels=levels)
    }
    x.excluded <- list()
    for (name in excluded.vars){
        levels <- mod$xlevels[[name]]
        fac <- !is.null(levels)
        level <- if (fac) levels[1] else mean(X[, name])
        if (fac) factor.levels[[name]] <- levels
        x.excluded[[name]] <- list(name=name, is.factor=fac,
                                   level=level)
    }
    dims <- sapply(x, function(x) length(x$levels))
    len <- prod(dims)
    n.basic <- length(basic.vars)
    n.excluded <- length(excluded.vars)
    n.vars <- n.basic + n.excluded
    predict.data <-matrix('', len, n.vars)
    excluded <- sapply(x.excluded, function(x) x$level)
    for (i in 1:len){
        subs <- subscripts(i, dims)
        for (j in 1:n.basic){
            predict.data[i,j] <- x[[j]]$levels[subs[j]]
        }
        if (n.excluded > 0)
            predict.data[i, (n.basic+1):n.vars] <- excluded
    }
    colnames(predict.data) <- c(sapply(x, function(x) x$name),
                                sapply(x.excluded, function(x) x$name))
    colclasses <- lapply(X, class)
    colclasses[colclasses == "matrix"] <- "numeric"
    predict.data <-  matrix.to.df(predict.data, colclasses=colclasses)
    list(predict.data=predict.data, factor.levels=factor.levels, 
         factor.cols=factor.cols, mod.aug=mod.aug, term=term, n.basic=n.basic,
         x=x, X.mod=X.mod, cnames=cnames, X=X)   
}

#fixup.model.matrix <- function(mod, mod.matrix, mod.matrix.all, X.mod, mod.aug, 
#	factor.cols, cnames, term, typical, given.values){
#	attr(mod.matrix, "assign") <- attr(mod.matrix.all, "assign")
#	stranger.cols <- factor.cols & 
#		apply(outer(strangers(term, mod, mod.aug), attr(mod.matrix,'assign'), '=='), 2, any)
#	if (has.intercept(mod)) stranger.cols[1] <- TRUE
#	if (any(stranger.cols)) {
#		mod.matrix[,stranger.cols] <- 
#			matrix(apply(as.matrix(X.mod[,stranger.cols]), 2, typical), 
#				nrow=nrow(mod.matrix), ncol=sum(stranger.cols), byrow=TRUE)
#		if (!is.null(given.values)){
#			stranger.names <- names(stranger.cols[stranger.cols])
#			given <- stranger.names %in% names(given.values)
#			if (any(given)) mod.matrix[,stranger.names[given]] <- given.values[stranger.names[given]]
#		} 
#	}
#	for (name in cnames){
#		components <- unlist(strsplit(name, ':'))
#		if (length(components) > 1) 
#			mod.matrix[,name] <- apply(mod.matrix[,components], 1, prod)
#	}
#	mod.matrix
#}

fixup.model.matrix <- function(mod, mod.matrix, mod.matrix.all, X.mod, mod.aug, 
	factor.cols, cnames, term, typical, given.values){
	attr(mod.matrix, "assign") <- attr(mod.matrix.all, "assign")
	stranger.cols <-  
		apply(outer(strangers(term, mod, mod.aug), attr(mod.matrix,'assign'), '=='), 2, any)
	if (has.intercept(mod)) stranger.cols[1] <- TRUE
	if (any(stranger.cols)) {
		facs <- factor.cols & stranger.cols
		covs <- (!factor.cols) & stranger.cols
		if (any(facs)) mod.matrix[,facs] <- 
				matrix(apply(as.matrix(X.mod[,facs]), 2, mean), 
					nrow=nrow(mod.matrix), ncol=sum(facs), byrow=TRUE)
		if (any(covs)) mod.matrix[,covs] <- 
				matrix(apply(as.matrix(X.mod[,covs]), 2, typical), 
					nrow=nrow(mod.matrix), ncol=sum(covs), byrow=TRUE)
		if (!is.null(given.values)){
			stranger.names <- cnames[stranger.cols]
			given <- stranger.names %in% names(given.values)
			if (any(given)) mod.matrix[,stranger.names[given]] <- 
					matrix(given.values[stranger.names[given]], nrow=nrow(mod.matrix), 
						ncol=length(stranger.names[given]), byrow=TRUE)
		} 
	}
	for (name in cnames){
		components <- unlist(strsplit(name, ':'))
		if (length(components) > 1) 
			mod.matrix[,name] <- apply(mod.matrix[,components], 1, prod)
	}
	mod.matrix
}

# the following function is a modification of code contributed by Steve Taylor

as.data.frame.eff <- function(x, row.names=NULL, optional=TRUE, transform=x$transformation$inverse, ...){
	result <- if (is.null(x$se)) data.frame(x$x, fit=transform(x$fit))
	else data.frame(x$x, fit=transform(x$fit), se=x$se, lower=transform(x$lower), upper=transform(x$upper))
    attr(result, "transformation") <- transform
    result
}

as.data.frame.effpoly <- function(x, row.names=NULL, optional=TRUE, ...){
	factors <- sapply(x$variables, function(x) x$is.factor)
	factor.levels <- lapply(x$variables[factors], function(x) x$levels)
	if (!length(factor.levels) == 0){
		factor.names <- names(factor.levels)
		for (fac in factor.names){
			x$x[[fac]] <- factor(x$x[[fac]], levels=factor.levels[[fac]])
		}
	}
	result <- data.frame(x$x, x$prob, x$logit)
	if (!is.null(x$confidence.level)) result <- cbind(result,
			x$se.prob, x$se.logit, x$lower.prob, x$upper.prob, x$lower.logit, x$upper.logit)
	result
}

as.data.frame.efflatent <- function(x, row.names=NULL, optional=TRUE, ...){
	if (is.null(x$se)) data.frame(x$x, fit=x$fit)
	else data.frame(x$x, fit=x$fit, se=x$se, lower=x$lower, upper=x$upper)
}

logit2p <- function(logit) 1/(1 + exp(-logit))

p2logit <- function(p) log(p/(1 - p))


lrug <- function(x) {
	if (length(unique(x)) < 0.8 * length(x)) x <- jitter(x)
	grid.segments(x, unit(0, "npc"), x, unit(0.5, "lines"),
		default.units="native")
}

#model.matrix.gls <- function(object, ...){
#	model.matrix(as.formula(object$call$model), data=eval(object$call$data))
##	model.matrix(lm(as.formula(object$call$model), data=eval(object$call$data)))
#}
#
#model.frame.gls <- function(formula, ...){
#	model.frame(as.formula(formula$call$model), data=eval(formula$call$data))
##	model.matrix(lm(as.formula(object$call$model), data=eval(object$call$data)))
#}
#
## model.response not generic
model.response.gls <- function(model){
	model.response(model.frame(as.formula(model$call$model), data=eval(model$call$data)))
}

## vcov method for eff objects

vcov.eff <- function(object, ...) object$vcov

## [ method for efflist objects

`[.efflist` <- function(x, ...){
	y <- NextMethod("[")
	class(y) <- class(x)
	y
}

### the following functions are for use by Effect() methods

Analyze.model <- function(focal.predictors, mod, xlevels, default.levels, formula.rhs){
    if ((!is.null(mod$na.action)) && class(mod$na.action) == "exclude") 
        class(mod$na.action) <- "omit"
    all.predictors <- all.vars(formula.rhs)
    check.vars <- !(focal.predictors %in% all.predictors)
    excluded.predictors <- setdiff(all.predictors, focal.predictors)
    number.bad <- sum(check.vars)
    if (any(check.vars)) {
        message <- if (number.bad == 1) paste("the following predictor is not in the model:", 
                                              all.predictors[check.vars])
        else paste("the following predictors are not in the model:", 
                   paste(all.predictors[check.vars], collapse=", "))
        stop(message)
    }
    X.mod <- model.matrix(mod)
    cnames <- colnames(X.mod)
    factor.cols <- rep(FALSE, length(cnames))
    names(factor.cols) <- cnames
    X <- model.frame(mod)
    for (name in all.predictors){
        if (is.factor(X[[name]])) factor.cols[grep(paste("^", name, sep=""), cnames)] <- TRUE
    }
    factor.cols[grep(":", cnames)] <- FALSE   
    X <- na.omit(expand.model.frame(mod, all.predictors))
    x <- list()
    factor.levels <- list()
    for (name in focal.predictors){
        levels <- mod$xlevels[[name]]
        fac <- !is.null(levels)
        if (!fac) {
            levels <- if (is.null(xlevels[[name]]))
                seq(min(X[, name]), max(X[,name]), length=default.levels)
            else xlevels[[name]]
        }
        else factor.levels[[name]] <- levels
        x[[name]] <- list(name=name, is.factor=fac, levels=levels)
    }
    x.excluded <- list()
    for (name in excluded.predictors){
        levels <- mod$xlevels[[name]]
        fac <- !is.null(levels)
        level <- if (fac) levels[1] else mean(X[, name])
        if (fac) factor.levels[[name]] <- levels
        x.excluded[[name]] <- list(name=name, is.factor=fac,
                                   level=level)
    }
    dims <- sapply(x, function(x) length(x$levels))
    len <- prod(dims)
    n.focal <- length(focal.predictors)
    n.excluded <- length(excluded.predictors)
    n.vars <- n.focal + n.excluded
    predict.data <-matrix('', len, n.vars)
    excluded <- sapply(x.excluded, function(x) x$level)
    for (i in 1:len){
        subs <- subscripts(i, dims)
        for (j in 1:n.focal){
            predict.data[i,j] <- x[[j]]$levels[subs[j]]
        }
        if (n.excluded > 0)
            predict.data[i, (n.focal + 1):n.vars] <- excluded
    }
    colnames(predict.data) <- c(sapply(x, function(x) x$name),
                                sapply(x.excluded, function(x) x$name))
    colclasses <- lapply(X, class)
    colclasses[colclasses == "matrix"] <- "numeric"
    predict.data <-  matrix.to.df(predict.data, colclasses=colclasses)
    list(predict.data=predict.data, factor.levels=factor.levels, 
         factor.cols=factor.cols, focal.predictors=focal.predictors, n.focal=n.focal,
         excluded.predictors=excluded.predictors, n.excluded=n.excluded,
         x=x, X.mod=X.mod, cnames=cnames, X=X)   
}

Fixup.model.matrix <- function(mod, mod.matrix, mod.matrix.all, X.mod,
		factor.cols, cnames, focal.predictors, excluded.predictors, typical, given.values){
	vars <- as.character(attr(terms(mod), "variables"))[-(1:2)]
	attr(mod.matrix, "assign") <- attr(mod.matrix.all, "assign")
	if (length(excluded.predictors) > 0){
		sel <- apply(sapply(excluded.predictors, matchVarName, expressions=vars), 1, any)
		strangers <- Strangers(mod, focal.predictors, excluded.predictors)
		stranger.cols <-  
				apply(outer(strangers, attr(mod.matrix,'assign'), '=='), 2, any)
	}
	else stranger.cols <- rep(FALSE, ncol(mod.matrix))
	if (has.intercept(mod)) stranger.cols[1] <- TRUE
	if (any(stranger.cols)) {
		facs <- factor.cols & stranger.cols
		covs <- (!factor.cols) & stranger.cols
		if (any(facs)) mod.matrix[,facs] <- 
					matrix(apply(as.matrix(X.mod[,facs]), 2, mean), 
							nrow=nrow(mod.matrix), ncol=sum(facs), byrow=TRUE)
		if (any(covs)) mod.matrix[,covs] <- 
					matrix(apply(as.matrix(X.mod[,covs]), 2, typical), 
							nrow=nrow(mod.matrix), ncol=sum(covs), byrow=TRUE)
		if (!is.null(given.values)){
			stranger.names <- cnames[stranger.cols]
			given <- stranger.names %in% names(given.values)
			if (any(given)) mod.matrix[,stranger.names[given]] <- 
						matrix(given.values[stranger.names[given]], nrow=nrow(mod.matrix), 
								ncol=length(stranger.names[given]), byrow=TRUE)
		} 
	}
	for (name in cnames){
		components <- unlist(strsplit(name, ':'))
		components <- components[components %in% cnames]
		if (length(components) > 1) 
			mod.matrix[,name] <- apply(mod.matrix[,components], 1, prod)
	}
	mod.matrix
}

matchVarName <- function(name, expressions){
	a <- !grepl(paste("[.]+", name, sep=""), expressions)
	b <- !grepl(paste(name, "[.]+", sep=""), expressions)
	c <- grepl(paste("\\b", name, "\\b", sep=""), expressions)
	a & b & c
}

Strangers <- function(mod, focal.predictors, excluded.predictors){
	names <- term.names(mod)
	if (has.intercept(mod)) names <- names[-1]
	sel <- apply(sapply(excluded.predictors, matchVarName, expressions=names), 1, any)
	(1:length(sel))[sel]
}

# the following is used by effect.multinom() and Effect.multinom()

eff.mul <- function(x0, B, se, m, p, r, V){
    mu <- exp(x0 %*% B)
    mu <- mu/(1 + sum(mu))
    mu[m] <- 1 - sum(mu)
    logits <- log(mu/(1 - mu))
    if (!se) return(list(p=mu, logits=logits))
    d <- array(0, c(m, m - 1, p))
    exp.x0.B <- as.vector(exp(x0 %*% B))
    sum.exp.x0.B <- sum(exp.x0.B)
    for (j in 1:(m-1)){
        d[m, j,] <- - exp.x0.B[j]*x0
        for (jj in 1:(m-1)){
            d[j, jj,] <- if (jj != j)
                - exp(x0 %*% (B[,jj] + B[,j]))*x0
            else exp.x0.B[j]*(1 + sum.exp.x0.B - exp.x0.B[j])*x0
        }
    }
    d <- d/(1 + sum.exp.x0.B)^2
    V.mu <- rep(0, m)
    for (j in 1:m){
        dd <- as.vector(t(d[j,,]))
        for (s in 1:r){
            for (t in 1:r){
                V.mu[j] <- V.mu[j] + V[s,t]*dd[s]*dd[t]
            }
        }
    }
    V.logits <- V.mu/(mu^2 * (1 - mu)^2)
    list(p=mu, std.err.p=sqrt(V.mu), logits=logits,
         std.error.logits=sqrt(V.logits))
}

# the following are used by effect.polr() and Effect.polr()

eff.polr <- function(x0, b, alpha, V, m, r, se){
    eta0 <- x0 %*% b
    mu <- rep(0, m)
    mu[1] <- 1/(1 + exp(alpha[1] + eta0))
    for (j in 2:(m-1)){
        mu[j] <- exp(eta0)*(exp(alpha[j - 1]) - exp(alpha[j]))/
            ((1 + exp(alpha[j - 1] + eta0))*(1 + exp(alpha[j] + eta0)))
    }
    mu[m] <- 1 - sum(mu)
    logits <- log(mu/(1 - mu))
    if (!se) return(list(p=mu, logits=logits))
    d <- matrix(0, m, r)
    d[1, 1] <- - exp(alpha[1] + eta0)/(1 + exp(alpha[1] + eta0))^2
    d[1, m:r] <- - exp(alpha[1] + eta0)*x0/(1 + exp(alpha[1] + eta0))^2
    for (j in 2:(m-1)){
        d[j, j-1] <- exp(alpha[j-1] + eta0)/(1 + exp(alpha[j-1] + eta0))^2
        d[j, j]   <- - exp(alpha[j] + eta0)/(1 + exp(alpha[j] + eta0))^2
        d[j, m:r] <- exp(eta0)*(exp(alpha[j]) - exp(alpha[j-1]))*
            (exp(alpha[j-1] + alpha[j] + 2*eta0) - 1) * x0 /
            (((1 + exp(alpha[j-1] + eta0))^2)*
            ((1 + exp(alpha[j] + eta0))^2))
    }
    d[m, m-1] <- exp(alpha[m-1] + eta0)/(1 + exp(alpha[m-1] + eta0))^2
    d[m, m:r] <- exp(alpha[m-1] + eta0)*x0/(1 + exp(alpha[m-1] + eta0))^2
    V.mu <- rep(0, m)
    for (j in 1:m){
        dd <- d[j,]
        for (s in 1:r){
            for (t in 1:r){
                V.mu[j] <- V.mu[j] + V[s,t]*dd[s]*dd[t]
            }
        }
    }
    V.logits <- V.mu/(mu^2 * (1 - mu)^2)
    list(p=mu, std.err.p=sqrt(V.mu), logits=logits,
         std.error.logits=sqrt(V.logits))
}

eff.latent <- function(X0, b, V, se){
    eta <- X0 %*% b
    if (!se) return(list(fit=eta))
    var <- diag(X0 %*% V %*% t(X0))
    list(fit=eta, se=sqrt(var))
}
