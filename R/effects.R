# effect generic and methods; allEffects
# John Fox and Jangman Hong
#  last modified 6 Februrary 2011 by J. Fox

effect <- function(term, mod, ...){
	UseMethod("effect", mod)
}

effect.lm <- function (term, mod, xlevels=list(), default.levels=10, given.values,
	se=TRUE, confidence.level=.95, 
	transformation=list(link=family(mod)$linkfun, inverse=family(mod)$linkinv), 
	typical=mean, ...){	
	if (missing(given.values)) given.values <- NULL
	else if (!all(which <- names(given.values) %in% names(coef(mod)))) 
		stop("given.values (", names(given.values[!which]),") not in the model")
	model.components <- analyze.model(term, mod, xlevels, default.levels)
	predict.data <- model.components$predict.data
	factor.levels <- model.components$factor.levels
	factor.cols <- model.components$factor.cols
	mod.aug <- model.components$mod.aug
	term <- model.components$term
	n.basic <- model.components$n.basic
	x <- model.components$x
	X.mod <- model.components$X.mod
	cnames <- model.components$cnames
	X <- model.components$X	
	formula.rhs <- formula(mod)[c(1,3)]  
	nrow.X <- nrow(X)
	mf <- model.frame(formula.rhs, data=rbind(X[,names(predict.data),drop=FALSE], predict.data), 
		xlev=factor.levels)
	mod.matrix.all <- model.matrix(formula.rhs, data=mf, contrasts.arg=mod$contrasts)
	mod.matrix <- mod.matrix.all[-(1:nrow.X),]
	fit.1 <- na.omit(predict(mod))
	wts <- mod$weights
	if (is.null(wts)) wts <- rep(1, length(fit.1))
	mod.2 <- lm.wfit(mod.matrix.all[1:nrow.X,], fit.1, wts)
	class(mod.2) <- "lm"
	y <- if(inherits(mod, "glm")) mod$y else na.omit(model.response(model.frame(mod)))
	discrepancy <- 100*mean(abs(fitted(mod.2)- fit.1)/(1e-10 + mean(abs(fit.1))))
	if (discrepancy > 1e-3) warning(paste("There is a discrepancy of", round(discrepancy, 3),
				"percent \n     in the 'safe' predictions used to generate effect", term))
	mod.matrix <- fixup.model.matrix(mod, mod.matrix, mod.matrix.all, X.mod, mod.aug, 
		factor.cols, cnames, term, typical, given.values)	
	effect <- mod.matrix %*% mod.2$coefficients
	result <- list(term=term, formula=formula(mod), response=response.name(mod),
		variables=x, fit=effect, 
		x=predict.data[,1:n.basic, drop=FALSE], model.matrix=mod.matrix, 
		data=X, discrepancy=discrepancy)
	if (se){
		if (any(family(mod)$family == c('binomial', 'poisson'))){
			dispersion <-  1
			z <- qnorm(1 - (1 - confidence.level)/2)
		}
		else {
			dispersion <- sum(wts * mod$residuals^2)/mod$df.residual
			z <- qt(1 - (1 - confidence.level)/2, df=mod$df.residual)
		}
		mod.2$terms <- mod$terms
		V <- dispersion * summary.lm(mod.2)$cov
		vcov <- mod.matrix %*% V %*% t(mod.matrix)
		rownames(vcov) <- colnames(vcov) <- NULL
		var <- diag(vcov)
		result$vcov <- vcov
		result$se <- sqrt(var)        
		result$lower <- effect - z*result$se
		result$upper <- effect + z*result$se
		result$confidence.level <- confidence.level
	}
	if (is.null(transformation$link) && is.null(transformation$inverse)){
		transformation$link <- I
		transformation$inverse <- I
	}
	result$transformation <- transformation
	class(result)<-'eff'
	result
}


effect.gls <- function (term, mod, xlevels=list(), default.levels=10, given.values,
		se=TRUE, confidence.level=.95, 
		transformation=NULL, 
		typical=mean, ...){	
	if (missing(given.values)) given.values <- NULL
	else if (!all(which <- names(given.values) %in% names(coef(mod)))) 
		stop("given.values (", names(given.values[!which]),") not in the model")
	mod.lm <- lm(as.formula(mod$call$model), data=eval(mod$call$data))
	model.components <- analyze.model(term, mod.lm, xlevels, default.levels)
	predict.data <- model.components$predict.data
	factor.levels <- model.components$factor.levels
	factor.cols <- model.components$factor.cols
	mod.aug <- model.components$mod.aug
	term <- model.components$term
	n.basic <- model.components$n.basic
	x <- model.components$x
	X.mod <- model.components$X.mod
	cnames <- model.components$cnames
	X <- model.components$X	
	formula.rhs <- formula(mod)[c(1,3)]  
	nrow.X <- nrow(X)
	mf <- model.frame(formula.rhs, data=rbind(X[,names(predict.data),drop=FALSE], predict.data), 
			xlev=factor.levels)
	mod.matrix.all <- model.matrix(formula.rhs, data=mf, contrasts.arg=mod$contrasts)
	mod.matrix <- mod.matrix.all[-(1:nrow.X),]
	fit.1 <- na.omit(predict(mod))
	mod.2 <- lm.fit(mod.matrix.all[1:nrow.X,], fit.1)
	class(mod.2) <- "lm"
	assign(".y", na.omit(model.response.gls(mod)), envir=.GlobalEnv)
	assign(".X", na.omit(mod.matrix.all[1:nrow.X,]), envir=.GlobalEnv)
	mod.3 <- update(mod, .y ~ .X - 1)
	remove(".X", ".y", envir=.GlobalEnv)
	discrepancy <- 100*mean(abs(fitted(mod.2)- fit.1)/(1e-10 + mean(abs(fit.1))))
	if (discrepancy > 1e-3) warning(paste("There is a discrepancy of", round(discrepancy, 3),
						"percent \n     in the 'safe' predictions used to generate effect", term))
	mod.matrix <- fixup.model.matrix(mod.lm, mod.matrix, mod.matrix.all, X.mod, mod.aug, 
			factor.cols, cnames, term, typical, given.values)	
	effect <- mod.matrix %*% mod.2$coefficients
	result <- list(term=term, formula=formula(mod), response=response.name(mod),
			variables=x, fit=effect, 
			x=predict.data[,1:n.basic, drop=FALSE], model.matrix=mod.matrix, 
			data=X, discrepancy=discrepancy)
	if (se){
		df.res <- mod$dims[["N"]] - mod$dims[["p"]]
		z <- qt(1 - (1 - confidence.level)/2, df=df.res)
		mod.2$terms <- terms(mod)
		V <- vcov(mod.3)
		vcov <- mod.matrix %*% V %*% t(mod.matrix)
		rownames(vcov) <- colnames(vcov) <- NULL
		var <- diag(vcov)
		result$vcov <- vcov
		result$se <- sqrt(var)        
		result$lower <- effect - z*result$se
		result$upper <- effect + z*result$se
		result$confidence.level <- confidence.level
	}
	if (is.null(transformation$link) && is.null(transformation$inverse)){
		transformation$link <- I
		transformation$inverse <- I
	}
	result$transformation <- transformation
	class(result)<-'eff'
	result
}


effect.multinom <- function(term, mod, 
	confidence.level=.95, xlevels=list(), default.levels=10, 
	given.values, se=TRUE, typical=mean, ...){	
	eff.mul <- function(x0){
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
	if (length(mod$lev) < 3) stop("effects for multinomial logit model only available for response levels > 2")
	if (missing(given.values)) given.values <- NULL
	else if (!all(which <- colnames(given.values) %in% names(coef(mod)))) 
		stop("given.values (", colnames(given.values[!which]),") not in the model")
	# refit model to produce 'safe' predictions when the model matrix includes
	#   terms -- e.g., poly(), bs() -- whose basis depends upon the data
	fit1 <- predict(mod, type="probs")
	model.components <- analyze.model(term, mod, xlevels, default.levels)
	predict.data <-model.components$predict.data
	factor.levels <- model.components$factor.levels
	factor.cols <- model.components$factor.cols
	mod.aug <- model.components$mod.aug
	term <- model.components$term
	n.basic <- model.components$n.basic
	x<- model.components$x
	X.mod <- model.components$X.mod
	cnames<- model.components$cnames
	X <- model.components$X	
	formula.rhs <- formula(mod)[c(1,3)]
	newdata <- predict.data
	newdata[[as.character(formula(mod)[2])]] <- rep(mod$lev[1], nrow(newdata))
	extras <- setdiff(all.vars(formula(mod)), names(model.frame(mod)))
	X <- if (length(extras) == 0) model.frame(mod)
		else {
			if (is.null(mod$call$data))
				mod$call$data <- environment(formula(mod))
			expand.model.frame(mod, extras)
		}
	X <- na.omit(X)
	nrow.X <- nrow(X)
	data <- rbind(X[,names(newdata),drop=FALSE], newdata)
	data$wt <- rep(0, nrow(data))
	data$wt[1:nrow.X] <- weights(mod)
	mod.matrix.all <- model.matrix(formula.rhs, data=data, contrasts.arg=mod$contrasts)
	X0 <- mod.matrix.all[-(1:nrow.X),]
	X0 <- fixup.model.matrix(mod, X0, mod.matrix.all, X.mod, mod.aug, factor.cols, 
		cnames, term, typical, given.values)
	resp.names <- make.names(mod$lev, unique=TRUE)
	resp.names <- c(resp.names[-1], resp.names[1]) # make the last level the reference level
	mod <- multinom(formula(mod), data=data, Hess=TRUE, weights=wt)	
	fit2 <- predict(mod, type="probs")[1:nrow.X,]
	fit1 <- na.omit(as.vector(p2logit(fit1)))
	fit2 <- as.vector(p2logit(fit2))
#	discrepancy <- 100*sqrt(mean((fit1 - fit2)^2)/mean(fit1^2))
	discrepancy <- 100*mean(abs(fit1 - fit2)/(1e-10 + mean(abs(fit1))))
	if (discrepancy > 0.1) warning(paste("There is a discrepancy of", round(discrepancy, 3),
				"percent \n     in the 'safe' predictions used to generate effect", term))
	B <- t(coef(mod))
	V <- vcov(mod)
	m <- ncol(B) + 1
	p <- nrow(B)
	r <- p*(m - 1)	
	n <- nrow(X0)
	P <- Logit <- matrix(0, n, m)
	colnames(P) <-  paste("prob.", resp.names, sep="")
	colnames(Logit) <-  paste("logit.", resp.names, sep="")
	if (se){
		z <- qnorm(1 - (1 - confidence.level)/2)
		Lower.P <- Upper.P <- Lower.logit <- Upper.logit <- SE.P <- SE.logit <- matrix(0, n, m)
		colnames(Lower.logit) <-  paste("L.logit.", resp.names, sep="")
		colnames(Upper.logit) <-  paste("U.logit.", resp.names, sep="")
		colnames(Lower.P) <-  paste("L.prob.", resp.names, sep="")
		colnames(Upper.P) <-  paste("U.prob.", resp.names, sep="")
		colnames(SE.P) <-  paste("se.prob.", resp.names, sep="")
		colnames(SE.logit) <-  paste("se.logit.", resp.names, sep="")
	}
	for (i in 1:n){
		res <- eff.mul(X0[i,]) # compute effects
		P[i,] <- prob <- res$p # fitted probabilities
		Logit[i,] <- logit <- res$logits # fitted logits
		if (se){
			SE.P[i,] <- se.p <- res$std.err.p # std. errors of fitted probs		
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
	result <- list(term=term, formula=formula(mod), response=response.name(mod),
		y.levels=mod$lev, variables=x, x=predict.data[,1:n.basic, drop=FALSE],
		model.matrix=X0, data=X, discrepancy=discrepancy, model="multinom",
		prob=P, logit=Logit)
	if (se) result <- c(result, list(se.prob=SE.P, se.logit=SE.logit,
				lower.logit=Lower.logit, upper.logit=Upper.logit, 
				lower.prob=Lower.P, upper.prob=Upper.P,
				confidence.level=confidence.level))
	class(result) <-'effpoly'
	result
}

effect.polr <- function(term, mod, 
	confidence.level=.95, xlevels=list(), default.levels=10, 
	given.values, se=TRUE, typical=mean, latent=FALSE, ...){
	if (mod$method != "logistic") stop('method argument to polr must be "logistic"')	
	if (missing(given.values)) given.values <- NULL
	else if (!all(which <- names(given.values) %in% names(coef(mod)))) 
		stop("given.values (", names(given.values[!which]),") not in the model")
	eff.polr <- function(x0){
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
	eff.latent <- function(X0, b, V){
		eta <- X0 %*% b
		if (!se) return(list(fit=eta))
		var <- diag(X0 %*% V %*% t(X0))
		list(fit=eta, se=sqrt(var))
	}
	# refit model to produce 'safe' predictions when the model matrix includes
	#   terms -- e.g., poly(), bs() -- whose basis depends upon the data
	fit1 <- predict(mod, type="probs")
	model.components <- analyze.model(term, mod, xlevels, default.levels)
	predict.data <- model.components$predict.data
	factor.levels <- model.components$factor.levels
	factor.cols <- model.components$factor.cols
	mod.aug <- model.components$mod.aug
	term <- model.components$term
	n.basic <- model.components$n.basic
	x <- model.components$x
	X.mod <- model.components$X.mod
	cnames <- model.components$cnames
#	X <- na.omit(model.components$X)
	formula.rhs <- formula(mod)[c(1,3)]
	newdata <- predict.data
	newdata[[as.character(formula(mod)[2])]] <- rep(mod$lev[1], nrow(newdata))
	extras <- setdiff(all.vars(formula(mod)), names(model.frame(mod)))
	X <- if (length(extras) == 0) model.frame(mod)
		else {
			if (is.null(mod$call$data))
				mod$call$data <- environment(formula(mod))
			expand.model.frame(mod, extras)
		}
	X <- na.omit(X)
	nrow.X <- nrow(X)
	data <- rbind(X[,names(newdata),drop=FALSE], newdata)
	wts <- mod$model[["(weights)"]]
	if (is.null(wts)) wts <- 1
	data$wt <- rep(0, nrow(data))
	data$wt[1:nrow.X] <- wts
	mod.matrix.all <- model.matrix(formula.rhs, data=data, contrasts.arg=mod$contrasts)
	X0 <- mod.matrix.all[-(1:nrow.X),]
	X0 <- fixup.model.matrix(mod, X0, mod.matrix.all, X.mod, mod.aug, factor.cols, 
		cnames, term, typical, given.values)
	resp.names <- make.names(mod$lev, unique=TRUE)
	mod <- polr(formula(mod), data=data, Hess=TRUE, weights=wt)
	fit2 <- predict(mod, type="probs")[1:nrow.X,]
	fit1 <- na.omit(as.vector(p2logit(fit1)))
	fit2 <- na.omit(as.vector(p2logit(fit2)))
#	discrepancy <- 100*sqrt(mean((fit1 - fit2)^2)/mean(fit1^2))
	discrepancy <- 100*mean(abs(fit1 - fit2)/(1e-10 + mean(abs(fit1))))
	if (discrepancy > 0.1) warning(paste("There is a discrepancy of", round(discrepancy, 3),
				"percent \n     in the 'safe' predictions used to generate effect", term))
	X0 <- X0[,-1, drop=FALSE]
	b <- coef(mod)
	p <- length(b)  # corresponds to p - 1 in the text
	alpha <- - mod$zeta  # intercepts are negatives of thresholds
	z <- qnorm(1 - (1 - confidence.level)/2)
	result <- list(term=term, formula=formula(mod), response=response.name(mod),
		y.levels=mod$lev, variables=x, 
		x=predict.data[,1:n.basic, drop=FALSE],
		model.matrix=X0, data=X, discrepancy=discrepancy, model="polr")
	if (latent){
		res <- eff.latent(X0, b, vcov(mod)[1:p, 1:p])
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
	V <- vcov(mod)[indices, indices]
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
		res <- eff.polr(X0[i,]) # compute effects
		P[i,] <- prob <- res$p # fitted probabilities
		Logit[i,] <- logit <- res$logits # fitted logits
		if (se){
			SE.P[i,] <- se.p <- res$std.err.p # std. errors of fitted probs		
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

allEffects <- function(mod, ...) UseMethod("allEffects")

allEffects.default <- function(mod, ...){
	high.order.terms <- function(mod){
		names <- term.names(mod)
		if (has.intercept(mod)) names<-names[-1]
		rel <- lapply(names, descendants, mod=mod)
		(1:length(names))[sapply(rel, function(x) length(x)==0)]
	}
	names <- term.names(mod)
	if (has.intercept(mod)) names <- names[-1]
	if (length(names) == 0) stop("the model contains no terms (beyond a constant)")
	terms <- names[high.order.terms(mod)]
	result <- lapply(terms, effect, mod=mod, ...)
	names(result) <- terms
	class(result) <- 'efflist'
	result
}

allEffects.gls <- function(mod, ...){
	high.order.terms <- function(mod){
		mod <- lm(as.formula(mod$call$model), data=eval(mod$call$data))
		names <- term.names(mod)
		if (has.intercept(mod)) names<-names[-1]
		rel <- lapply(names, descendants, mod=mod)
		(1:length(names))[sapply(rel, function(x) length(x)==0)]
	}
	names <- term.names(mod)
	if (has.intercept(mod)) names <- names[-1]
	if (length(names) == 0) stop("the model contains no terms (beyond a constant)")
	terms <- names[high.order.terms(mod)]
	result <- lapply(terms, effect, mod=mod, ...)
	names(result) <- terms
	class(result) <- 'efflist'
	result
}

all.effects <- function(...){
	.Deprecated("allEffects")
	allEffects(...)
}