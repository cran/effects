# plot.eff method for effects package, moved here from plot-summary-print-methods.R
# The plot.effpoly method remains there for now.

# the following functions aren't exported

make.ticks <- function(range, link, inverse, at, n) {
	warn <- options(warn=-1)
	on.exit(options(warn))
	link <- if (is.null(link)) 
				function(x) nlm(function(y) (inverse(y) - x)^2, 
							mean(range))$estimate
			else link
	if (is.null(n)) n <- 5
	labels <- if (is.null(at)){
				labels <- pretty(sapply(range, inverse), n=n+1)
			}
			else at
	ticks <- sapply(labels, link)
	list(at=ticks, labels=format(labels))
}

range.adj <- function(x){
	range <- range(x, na.rm=TRUE)
	c(range[1] - .025*(range[2] - range[1]),                                              
			range[2] + .025*(range[2] - range[1]))
}

# added, modified from http://www.r-bloggers.com/confidence-bands-with-lattice-and-r/

panel.bands <- function(x, y, upper, lower, fill, col,
		subscripts, ..., font, fontface)
{
	if(!missing(subscripts)) {
		upper <- upper[subscripts]
		lower <- lower[subscripts]
	}
	panel.polygon(c(x, rev(x)), c(upper, rev(lower)),
			col = fill, fill=fill, border = FALSE,
			...)
}


# modified by Michael Friendly: added key.args:
# modified by Michael Friendly: added ci.style="bands"
# modified by Michael Friendly: added lwd= argument for llines (not used elsewhere)
# modified by Michael Friendly: added alpha= argument for ci.style="bands"

plot.eff <- function(x, x.var=which.max(levels),
		z.var=which.min(levels), multiline=is.null(x$se), rug=TRUE, xlab,
		ylab, main=paste(effect, "effect plot"),
		colors=palette(), symbols=1:10, lines=1:10, cex=1.5, lwd=2, ylim, xlim=NULL,
		factor.names=TRUE, ci.style, alpha=0.3, 
		type=c("response", "link"), ticks=list(at=NULL, n=5),  
		alternating=TRUE, rotx=0, roty=0, grid=FALSE, layout, rescale.axis=TRUE, 
		transform.x=NULL, ticks.x=NULL,
		key.args=NULL, 
		row=1, col=1, nrow=1, ncol=1, more=FALSE, ...)
{  
	ci.style <- if(missing(ci.style)) NULL else 
				match.arg(ci.style, c("bars", "lines", "bands", "none")) 
	type <- match.arg(type)
	thresholds <- x$thresholds
	has.thresholds <- !is.null(thresholds)
	if (missing(ylab)){
		ylab <- if (has.thresholds) paste(x$response, ": ", paste(x$y.levels, collapse=", "), sep="")
				else x$response
	}     
	if (has.thresholds){ 
		threshold.labels <- abbreviate(x$y.levels, minlength=1)
		threshold.labels <- paste(" ", 
				paste(threshold.labels[-length(threshold.labels)], threshold.labels[-1], sep=" - "),
				" ", sep="")
	}
	trans.link <- x$transformation$link
	trans.inverse <- x$transformation$inverse
	if (!rescale.axis){
		x$lower[!is.na(x$lower)] <- trans.inverse(x$lower[!is.na(x$lower)])
		x$upper[!is.na(x$upper)] <- trans.inverse(x$upper[!is.na(x$upper)])
		x$fit[!is.na(x$fit)] <- trans.inverse(x$fit)[!is.na(x$fit)]
		trans.link <- trans.inverse <- I
	}
	# require(lattice)
	split <- c(col, row, ncol, nrow)
	ylab # force evaluation
	x.data <- x$data
	effect <- paste(sapply(x$variables, "[[", "name"), collapse="*")
	vars <- x$variables
	x <- as.data.frame(x, transform=I)
	for (i in 1:length(vars)){
		if (!(vars[[i]]$is.factor)) next
		x[,i] <- factor(x[,i], levels=vars[[i]]$levels)
	}
	has.se <- !is.null(x$se)
	n.predictors <- ncol(x) - 1 - 3*has.se
	if (n.predictors == 1){
		### factor no other predictors
		if (is.factor(x[,1])){
			ci.style <- if(is.null(ci.style)) "bars" else ci.style
			range <- if(has.se & ci.style!="none")
						range(c(x$lower, x$upper), na.rm=TRUE) else range(x$fit, na.rm=TRUE)
			ylim <- if (!missing(ylim)) ylim else c(range[1] - .025*(range[2] - range[1]),                                              
								range[2] + .025*(range[2] - range[1]))
			tickmarks <- if (type == "response") make.ticks(ylim, 
								link=trans.link, inverse=trans.inverse, at=ticks$at, n=ticks$n)
					else make.ticks(ylim, link=I, inverse=I, at=ticks$at, n=ticks$n)
			levs <- levels(x[,1])  
			plot <- xyplot(eval(parse(
									text=paste("fit ~ as.numeric(", names(x)[1], ")"))), 
					strip=function(...) strip.default(..., strip.names=c(factor.names, TRUE)),
					panel=function(x, y, lower, upper, has.se, ...){
						if (grid) panel.grid()
						good <- !is.na(y)
						if (has.se){ 
							if (ci.style == "bars"){
								larrows(x0=x[good], y0=lower[good], x1=x[good], y1=upper[good], angle=90, 
										code=3, col=colors[2], length=0.125*cex/1.5)
							}
							else if(ci.style == "lines") {
								llines(x[good], lower[good], lty=2, col=colors[2])
								llines(x[good], upper[good], lty=2, col=colors[2])
							}
							else{ if(ci.style == "bands") {
									panel.bands(x[good], y[good], upper[good], lower[good], fill=colors[2], alpha=alpha)
								}}
						}
						llines(x[good], y[good], lwd=lwd, col=colors[1], type='b', pch=19, cex=cex, ...)
						if (has.thresholds){
							panel.abline(h=thresholds, lty=3)
							panel.text(rep(current.panel.limits()$xlim[1], length(thresholds)), 
									thresholds, threshold.labels, adj=c(0,0), cex=0.75)
							panel.text(rep(current.panel.limits()$xlim[2], length(thresholds)), 
									thresholds, threshold.labels, adj=c(1,0), cex=0.75)
						}
					},
					ylim=ylim,
					ylab=ylab,
					xlab=if (missing(xlab)) names(x)[1] else xlab,
					scales=list(x=list(at=1:length(levs), labels=levs, rot=rotx), 
							y=list(at=tickmarks$at, labels=tickmarks$labels, rot=roty),
							alternating=alternating, y=roty),
					main=main,
					lower=x$lower, upper=x$upper, has.se=has.se, data=x, ...)
			result <- update(plot, layout = if (missing(layout)) c(0, prod(dim(plot))) 
							else layout)
			result$split <- split
			result$more <- more
			class(result) <- c("plot.eff", class(result))
		}  
		### variate, no other predictors      
		else {
			ci.style <- if(is.null(ci.style)) "lines" else ci.style
			range <- if(has.se & ci.style!="none")
						range(c(x$lower, x$upper), na.rm=TRUE) else range(x$fit, na.rm=TRUE)
			ylim <- if (!missing(ylim)) ylim else c(range[1] - .025*(range[2] - range[1]),                                              
								range[2] + .025*(range[2] - range[1]))
			tickmarks <- if (type == "response") make.ticks(ylim, 
								link=trans.link, inverse=trans.inverse, at=ticks$at, n=ticks$n)
					else make.ticks(ylim, link=I, inverse=I, at=ticks$at, n=ticks$n)
			nm <- names(x)[1]
			x.vals <- x.data[, nm]   
			if (nm %in% names(ticks.x)){
				at <- ticks.x[[nm]]$at
				n <- ticks.x[[nm]]$n
			}
			else{
				at <- NULL
				n <- 5
			}
			xlm <- if (nm %in% names(xlim)){
						xlim[[nm]]
					}
					else range.adj(x[nm]) # range(x.vals)
			tickmarks.x <- if ((nm %in% names(transform.x)) && !(is.null(transform.x))){
						trans <- transform.x[[nm]]$trans
#                make.ticks(trans(xlm), link=transform.x[[nm]]$trans, inverse=transform.x[[nm]]$inverse, at=ticks.x$at, n=ticks.x$n)
						make.ticks(trans(xlm), link=transform.x[[nm]]$trans, inverse=transform.x[[nm]]$inverse, at=at, n=n)
					}
					else {
						trans <- I
						make.ticks(xlm, link=I, inverse=I, at=at, n=n)
					}
			plot <- xyplot(eval(parse(
									text=paste("fit ~ trans(", names(x)[1], ")"))),
					strip=function(...) strip.default(..., strip.names=c(factor.names, TRUE)),
					panel=function(x, y, x.vals, rug, lower, upper, has.se, ...){
						if (grid) panel.grid()
						good <- !is.na(y)
						axis.length <- diff(range(x))
						llines(x[good], y[good], lwd=lwd, col=colors[1], ...)
						if (rug) lrug(trans(x.vals))
						if (has.se){  
							if (ci.style == "bars"){
								larrows(x0=x[good], y0=lower[good], 
										x1=x[good], y1=upper[good], 
										angle=90, code=3, col=eval(colors[2]), 
										length=.125*cex/1.5)
							}
							else if(ci.style == "lines") {
								llines(x[good], lower[good], lty=2, col=colors[2])
								llines(x[good], upper[good], lty=2, col=colors[2])
							}
							else{ if(ci.style == "bands") {
									panel.bands(x[good], y[good], upper[good], lower[good], fill=colors[2], alpha=alpha)
								}}
						}
						if (has.thresholds){
							panel.abline(h=thresholds, lty=3)
							panel.text(rep(current.panel.limits()$xlim[1], length(thresholds)), 
									thresholds, threshold.labels, adj=c(0,0), cex=0.75)
							panel.text(rep(current.panel.limits()$xlim[2], length(thresholds)), 
									thresholds, threshold.labels, adj=c(1,0), cex=0.75)
						}
					},
					ylim=ylim,
					xlim=suppressWarnings(trans(xlm)),
					ylab=ylab,
					xlab=if (missing(xlab)) names(x)[1] else xlab,
					x.vals=x.vals, rug=rug,
					main=main,
					lower=x$lower, upper=x$upper, has.se=has.se, data=x, 
					scales=list(y=list(at=tickmarks$at, labels=tickmarks$labels, rot=roty),
							x=list(at=tickmarks.x$at, labels=tickmarks.x$labels, rot=rotx), alternating=alternating), ...)
			result <- update(plot, layout = if (missing(layout)) c(0, prod(dim(plot))) 
							else layout)
			result$split <- split
			result$more <- more
			class(result) <- c("plot.eff", class(result))
		}
		return(result)
	}
	###  more than one variate
	predictors <- names(x)[1:n.predictors]
	levels <- sapply(apply(x[,predictors], 2, unique), length)
	if (is.character(x.var)) {
		which.x <- which(x.var == predictors)
		if (length(which.x) == 0) stop(paste("x.var = '", x.var, "' is not in the effect.", sep=""))
		x.var <- which.x
	}
	if (is.character(z.var)) {
		which.z <- which(z.var == predictors)
		if (length(which.z) == 0) stop(paste("z.var = '", z.var, "' is not in the effect.", sep=""))
		z.var <- which.z
	}    
	if (x.var == z.var) z.var <- z.var + 1
	### multiline
	if (multiline){ 
		ci.style <- if(is.null(ci.style)) "none" else ci.style
		if(ci.style == "lines") { 
			cat("Confidence interval style 'lines' changed to 'bars'\n")
			ci.style <- "bars"}
		range <- if (has.se && ci.style !="none")
					range(c(x$lower, x$upper), na.rm=TRUE) else range(x$fit, na.rm=TRUE)
		ylim <- if (!missing(ylim)) ylim else c(range[1] - .025*(range[2] - range[1]),                                              
							range[2] + .025*(range[2] - range[1]))
		tickmarks <- if (type == "response") make.ticks(ylim, link=trans.link, 
							inverse=trans.inverse, at=ticks$at, n=ticks$n)
				else make.ticks(ylim, link=I, inverse=I, at=ticks$at, n=ticks$n)
		zvals <- unique(x[, z.var])
		if (length(zvals) > min(c(length(colors), length(lines), length(symbols))))
			stop(paste('Not enough colors, lines, or symbols to plot', length(zvals), 'lines'))
		
		### multiline factor
		if (is.factor(x[,x.var])){
			levs <- levels(x[,x.var])
			key <- list(title=predictors[z.var], cex.title=1, border=TRUE,
					text=list(as.character(zvals)), 
					lines=list(col=colors[1:length(zvals)], lty=lines[1:length(zvals)], lwd=lwd), 
					points=list(col=colors[1:length(zvals)], pch=symbols[1:length(zvals)]))
			key <- c(key, key.args)
			plot <- xyplot(eval(parse( 
									text=paste("fit ~ as.numeric(", predictors[x.var], ")",
											if (n.predictors > 2) paste(" |", 
														paste(predictors[-c(x.var, z.var)])), collapse="*"))),
					strip=function(...) strip.default(..., strip.names=c(factor.names, TRUE)),
					panel=function(x, y, subscripts, z, lower, upper, show.se, ...){
						if (grid) panel.grid()
						for (i in 1:length(zvals)){
							sub <- z[subscripts] == zvals[i]
							good <- !is.na(y[sub])
							os <- if(show.se)
										(i - (length(zvals) + 1)/2) * (2/(length(zvals)-1)) * 
												.01 * (length(zvals) - 1) else 0
							llines(x[sub][good]+os, y[sub][good], lwd=lwd, type='b', col=colors[i], 
									pch=symbols[i], lty=lines[i], cex=cex, ...)
							if (show.se){
								larrows(x0=x[sub][good]+os, y0=lower[subscripts][sub][good], 
										x1=x[sub][good]+os, y1=upper[subscripts][sub][good], 
										angle=90, code=3, col=eval(colors[i]), 
										length=.125*cex/1.5)
							}
						} 
						if (has.thresholds){
							panel.abline(h=thresholds, lty=3)
							panel.text(rep(current.panel.limits()$xlim[1], length(thresholds)), 
									thresholds, threshold.labels, adj=c(0,0), cex=0.75)
							panel.text(rep(current.panel.limits()$xlim[2], length(thresholds)), 
									thresholds, threshold.labels, adj=c(1,0), cex=0.75)
						}
					},        
					ylim=ylim,
					ylab=ylab,
					xlab=if (missing(xlab)) predictors[x.var] else xlab,
					z=x[,z.var],
					scales=list(x=list(at=1:length(levs), labels=levs, rot=rotx), 
							y=list(at=tickmarks$at, labels=tickmarks$labels, rot=roty),
							alternating=alternating),
					zvals=zvals,
					main=main,
					key=key,
					lower=x$lower, upper=x$upper, 
					show.se=has.se && ci.style=="bars", 
					data=x, ...)
			result <- update(plot, layout = if (missing(layout)) 
								c(0, prod(dim(plot))) else layout)
			result$split <- split
			result$more <- more
			class(result) <- c("plot.eff", class(result))
		} 
		### multiline variate   
		else{
			nm <- names(x)[x.var]
			x.vals <- x.data[, nm]   
			if (nm %in% names(ticks.x)){
				at <- ticks.x[[nm]]$at
				n <- ticks.x[[nm]]$n
			}
			else{
				at <- NULL
				n <- 5
			}
			xlm <- if (nm %in% names(xlim)){
						xlim[[nm]]
					}
					else range.adj(x[nm]) # range(x.vals)
			tickmarks.x <- if ((nm %in% names(transform.x)) && !(is.null(transform.x))){
						trans <- transform.x[[nm]]$trans
						# make.ticks(trans(xlm), link=transform.x[[nm]]$trans, inverse=transform.x[[nm]]$inverse, at=ticks.x$at, n=ticks.x$n)
						make.ticks(trans(xlm), link=transform.x[[nm]]$trans, inverse=transform.x[[nm]]$inverse, at=at, n=n)
					}
					else {
						trans <- I
						make.ticks(xlm, link=I, inverse=I, at=at, n=n)
					}
			key<-list(title=predictors[z.var], cex.title=1, border=TRUE,
					text=list(as.character(zvals)), 
					lines=list(col=colors[1:length(zvals)], lty=lines[1:length(zvals)], lwd=lwd))
			key <- c(key, key.args) 
			plot <- xyplot(eval(parse( 
									text=paste("fit ~trans(", predictors[x.var], ")", 
											if (n.predictors > 2) paste(" |", 
														paste(predictors[-c(x.var, z.var)])), collapse="*"))),
					strip=function(...) strip.default(..., strip.names=c(factor.names, TRUE)),
					# old           panel=function(x, y, subscripts, x.vals, rug, z, ...){
					panel=function(x, y, subscripts, x.vals, rug, z, lower, upper, show.se, ...){
						if (grid) panel.grid()
						if (rug) lrug(trans(x.vals))
						#new
						axis.length <- diff(range(x))
						#end 
						for (i in 1:length(zvals)){
							sub <- z[subscripts] == zvals[i]
							good <- !is.na(y[sub])
							llines(x[sub][good], y[sub][good], lwd=lwd, type='l', col=colors[i], lty=lines[i], cex=cex, ...)
							if(show.se){
								os <- (i - (length(zvals) + 1)/2) * (2/(length(zvals)-1)) * 
										.01 * axis.length
								larrows(x0=x[sub][good]+os, y0=lower[subscripts][sub][good], 
										x1=x[sub][good]+os, y1=upper[subscripts][sub][good], 
										angle=90, code=3, col=eval(colors[i]), 
										length=.125*cex/1.5)
							}
						}
						if (has.thresholds){
							panel.abline(h=thresholds, lty=3)
							panel.text(rep(current.panel.limits()$xlim[1], length(thresholds)), 
									thresholds, threshold.labels, adj=c(0,0), cex=0.75)
							panel.text(rep(current.panel.limits()$xlim[2], length(thresholds)), 
									thresholds, threshold.labels, adj=c(1,0), cex=0.75)
						}
					},
					ylim=ylim,
					xlim=suppressWarnings(trans(xlm)), 
					ylab=ylab,
					xlab=if (missing(xlab)) predictors[x.var] else xlab,
					x.vals=x.vals, rug=rug,
					z=x[,z.var],
					zvals=zvals,
					main=main,
					key=key, 
					#
					lower=x$lower, upper=x$upper, 
					show.se=has.se && ci.style =="bars",
					#
					data=x, scales=list(y=list(at=tickmarks$at, labels=tickmarks$labels),
							rot=roty, x=list(at=tickmarks.x$at, labels=tickmarks.x$labels, rot=rotx), 
							alternating=alternating),  ...)
			result <- update(plot, layout = if (missing(layout)) c(0, prod(dim(plot))) 
							else layout)
			result$split <- split
			result$more <- more
			class(result) <- c("plot.eff", class(result))
		}
		return(result)
	} 
	# multiplot factor
	ci.style <- if(is.null(ci.style)){
				if(is.factor(x[, x.var])) "bars" else "lines"} else ci.style
	range <- if (has.se && ci.style !="none")
				range(c(x$lower, x$upper), na.rm=TRUE) else range(x$fit, na.rm=TRUE)
	ylim <- if (!missing(ylim)) ylim else c(range[1] - .025*(range[2] - range[1]),                                              
						range[2] + .025*(range[2] - range[1]))
	tickmarks <- if (type == "response") make.ticks(ylim, link=trans.link, 
						inverse=trans.inverse, at=ticks$at, n=ticks$n)
			else make.ticks(ylim, link=I, inverse=I, at=ticks$at, n=ticks$n)  
	if (is.factor(x[,x.var])){
		levs <- levels(x[,x.var])
		plot <- xyplot(eval(parse( 
								text=paste("fit ~ as.numeric(", predictors[x.var], ") |", 
										paste(predictors[-x.var], collapse="*")))),
				strip=function(...) strip.default(..., strip.names=c(factor.names, TRUE)),
				panel=function(x, y, subscripts, lower, upper, has.se, ...){  
					if (grid) panel.grid()
					good <- !is.na(y)
					if (has.se){
						if (ci.style == "bars"){
							larrows(x0=x[good], y0=lower[subscripts][good], x1=x[good], y1=upper[subscripts][good], 
									angle=90, code=3, col=colors[2], length=0.125*cex/1.5)
						}
						else if(ci.style == "lines") {
							llines(x[good], lower[subscripts][good], lty=2, col=colors[2])
							llines(x[good], upper[subscripts][good], lty=2, col=colors[2])
						}
						else{ if(ci.style == "bands") {
								panel.bands(x[good], y[good], upper[subscripts][good], lower[subscripts][good], fill=colors[2], alpha=alpha)
							}}
					}
					llines(x[good], y[good], lwd=lwd, type='b', col=colors[1], pch=19, cex=cex, ...)
					if (has.thresholds){
						panel.abline(h=thresholds, lty=3)
						panel.text(rep(current.panel.limits()$xlim[1], length(thresholds)), 
								thresholds, threshold.labels, adj=c(0,0), cex=0.75)
						panel.text(rep(current.panel.limits()$xlim[2], length(thresholds)), 
								thresholds, threshold.labels, adj=c(1,0), cex=0.75)
					}
				},
				ylim=ylim,
				ylab=ylab,
				xlab=if (missing(xlab)) predictors[x.var] else xlab,
				scales=list(x=list(at=1:length(levs), labels=levs, rot=rotx), 
						y=list(at=tickmarks$at, labels=tickmarks$labels, rot=roty),
						alternating=alternating),
				main=main,
				lower=x$lower, upper=x$upper, has.se=has.se, data=x, ...)
		result <- update(plot, layout = if (missing(layout)) c(0, prod(dim(plot))) else layout)
		result$split <- split
		result$more <- more
		class(result) <- c("plot.eff", class(result))
	} 
	### multiplot variate   
	else{
		nm <- names(x)[x.var]
		x.vals <- x.data[, nm]   
		if (nm %in% names(ticks.x)){
			at <- ticks.x[[nm]]$at
			n <- ticks.x[[nm]]$n
		}
		else{
			at <- NULL
			n <- 5
		}
		xlm <- if (nm %in% names(xlim)){
					xlim[[nm]]
				}
				else range.adj(x[nm]) # range(x.vals)
		tickmarks.x <- if ((nm %in% names(transform.x)) && !(is.null(transform.x))){
					trans <- transform.x[[nm]]$trans
					#  make.ticks(trans(xlm), link=transform.x[[nm]]$trans, inverse=transform.x[[nm]]$inverse, at=ticks.x$at, n=ticks.x$n)
					make.ticks(trans(xlm), link=transform.x[[nm]]$trans, inverse=transform.x[[nm]]$inverse, at=at, n=n)
				}
				else {
					trans <- I
					make.ticks(xlm, link=I, inverse=I, at=at, n=n)
				}
		plot <- xyplot(eval(parse( 
								text=paste("fit ~ trans(", predictors[x.var], ") |", 
										paste(predictors[-x.var], collapse="*")))),
				strip=function(...) strip.default(..., strip.names=c(factor.names, TRUE)),
				panel=function(x, y, subscripts, x.vals, rug, lower, upper, has.se, ...){
					if (grid) panel.grid()
					good <- !is.na(y)
					llines(x[good], y[good], lwd=lwd, col=colors[1], ...)
					if (rug) lrug(trans(x.vals))
					if (has.se){  
						if (ci.style == "bars"){
							larrows(x0=x[good], y0=lower[subscripts][good], 
									x1=x[good], y1=upper[subscripts][good], 
									angle=90, code=3, col=eval(colors[2]), 
									length=.125*cex/1.5)
						}
						else if(ci.style == "lines") {
							llines(x[good], lower[subscripts][good], lty=2, col=colors[2])
							llines(x[good], upper[subscripts][good], lty=2, col=colors[2])
						}
						else{ if(ci.style == "bands") {
								panel.bands(x[good], y[good], upper[subscripts][good], lower[subscripts][good], fill=colors[2], alpha=alpha)
							}}
					}
					if (has.thresholds){
						panel.abline(h=thresholds, lty=3)
						panel.text(rep(current.panel.limits()$xlim[1], length(thresholds)), 
								thresholds, threshold.labels, adj=c(0,0), cex=0.75)
						panel.text(rep(current.panel.limits()$xlim[2], length(thresholds)), 
								thresholds, threshold.labels, adj=c(1,0), cex=0.75)
					}
				},
				ylim=ylim,
				xlim=suppressWarnings(trans(xlm)),
				ylab=ylab,
				xlab=if (missing(xlab)) predictors[x.var] else xlab,
				x.vals=x.vals, rug=rug,
				main=main,
				lower=x$lower, upper=x$upper, has.se=has.se, data=x, 
				scales=list(y=list(at=tickmarks$at, labels=tickmarks$labels, rot=roty),
						x=list(at=tickmarks.x$at, labels=tickmarks.x$labels, rot=rotx), 
						alternating=alternating), ...)
		result <- update(plot, layout = if (missing(layout)) c(0, prod(dim(plot))) else layout)
		result$split <- split
		result$more <- more
		class(result) <- c("plot.eff", class(result))
	}
	return(result)
}

print.plot.eff <- function(x, ...){
	NextMethod(split=x$split, more=x$more, ...)
	invisible(x)
}

plot.efflist <- function(x, selection, rows, cols, ask=FALSE, graphics=TRUE, ...){
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
				print(plot(x[[(i-1)*cols + j]], row=i, col=j, nrow=rows, ncol=cols, more=more, ...))
			}
		}
	}
}

