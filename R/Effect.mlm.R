#  Calculate Effects for term(s) in a Multivariate Linear Model 


Effect.mlm <- function(focal.predictors, mod, response, ...) {
	if (missing(response)) {
		mod.frame <- model.frame(mod)
    response <- colnames(model.response(mod.frame))
	}
	else if (is.numeric(response)) {
		mod.frame <- model.frame(mod)
    response.names <- colnames(model.response(mod.frame))
    response <- response.names[response]
	}
	
	if (length(response)==1) {
			mod.1 <- update(mod, as.formula(paste(response, " ~ .")))
			result <- Effect(focal.predictors, mod.1,  ...)
	}
	else {
		result <- as.list(NULL)
		for (resp in response) {
			mod.1 <- update(mod, as.formula(paste(resp, " ~ .")))
			lab <- resp
			result[[lab]] <- Effect(focal.predictors, mod.1,  ...)
		}
		class(result) <- "efflist"
	}
	result
}

