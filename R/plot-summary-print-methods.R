# plot, summary, and print methods for effects package
# John Fox and Jangman Hong
#  last modified 2012-11-30 by J. Fox
#  29 June 2011 added grid, rotx and roty arguments to the two plot methods
#   by S. Weisberg
#  21 Dec 2012 modest modification of empty cells with crossed factors
#  2013-01-17: Added factor.ci.style arg to plot.eff() and plot.effpoly(). J. Fox
#  2013-01-18: Added CI bars to multiline plots with factor.ci.style="bars"
#  2013-01-19: Renamed 'factor.ci.style' to 'ci.style'.  Added a 'none' option
#   extended to variate terms if multiline=TRUE, ci.style="bars"
#  2013-01-30: scale arrow "heads" for error bars relative to cex


summary.eff <- function(object, type=c("response", "link"), ...){
	result <- list()
	result$header <- paste("\n", gsub(":", "*", object$term), 'effect\n')
	result$offset <- object$offset
	type <- match.arg(type)
	if (type == "response") {
		object$fit <- object$transformation$inverse(object$fit)
		if (!is.null(object$confidence.level)){
			object$lower <- object$transformation$inverse(object$lower)
			object$upper <- object$transformation$inverse(object$upper)
		}
	}
	result$effect <- array(object$fit,     
			dim=sapply(object$variables, function(x) length(x$levels)),
			dimnames=lapply(object$variables, function(x) x$levels))
	if (!is.null(object$se)){
		result$lower.header <- paste('\n Lower', round(100*object$confidence.level, 2), 
				'Percent Confidence Limits\n')
		result$lower <- array(object$lower,   
				dim=sapply(object$variables, function(x) length(x$levels)),
				dimnames=lapply(object$variables, function(x) x$levels))
		result$upper.header <- paste('\n Upper', round(100*object$confidence.level, 2),
				'Percent Confidence Limits\n')
		result$upper <- array(object$upper,   
				dim=sapply(object$variables, function(x) length(x$levels)),
				dimnames=lapply(object$variables, function(x) x$levels))
	}
	if (object$discrepancy > 1e-3) result$warning <- paste("\nWarning: There is an average discrepancy of", 
				round(object$discrepancy, 3),
				"percent \n     in the 'safe' predictions for effect", object$term, '\n')
	class(result) <- "summary.eff"
	result
}

print.summary.eff <- function(x, ...){
	cat(x$header)
	if (x$offset != 0) cat("\noffset = ", x$offset, "\n\n")
	print(x$effect, ...)
	if (!is.null(x$lower)){
		cat(x$lower.header)
		print(x$lower, ...)
		cat(x$upper.header)
		print(x$upper, ...)
	}
	if (!is.null(x$thresholds)){
		cat("\nThresholds:\n")
		print(x$thresholds, ...)
	}
	if (!is.null(x$warning)) cat(x$warning)
	invisible(x)
}

print.eff <- function(x, type=c("response", "link"), ...){
	cat(paste("\n", gsub(":", "*", x$term), 'effect\n'))
	if (x$offset != 0) cat("\noffset = ", x$offset, "\n\n")
	type <- match.arg(type)
	if (type == "response") x$fit <- x$transformation$inverse(x$fit)
	table <- array(x$fit,     
			dim=sapply(x$variables, function(x) length(x$levels)),
			dimnames=lapply(x$variables, function(x) x$levels))
	print(table, ...)
	if (x$discrepancy > 1e-3) cat(paste("\nWarning: There is an average discrepancy of", 
						round(x$discrepancy, 3),
						"percent \n     in the 'safe' predictions for effect", x$term, '\n'))
	invisible(x)
}

print.efflist <- function(x, ...){
	cat(" model: ")
	form <- x[[1]]$formula
	attributes(form) <- NULL
	print(form)
	for (effect in names(x)){
		print(x[[effect]], ...)
	}
	invisible(x) 
}

summary.efflist <- function(object, ...){
	cat(" model: ")
	form <- object[[1]]$formula
	attributes(form) <- NULL
	print(form)
	for (effect in names(object)){
		print(summary(object[[effect]], ...))
	}
	invisible(NULL) 
}

# the following two functions aren't exported

make.ticks <- function(range, link, inverse, at, n) {
    warn <- options(warn=-1)
    on.exit(warn)
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

# modified by Michael Friendly: added key.args:

plot.eff <- function(x, x.var=which.max(levels),
    z.var=which.min(levels), multiline=is.null(x$se), rug=TRUE, xlab,
    ylab, main=paste(effect, "effect plot"),
    colors=palette(), symbols=1:10, lines=1:10, cex=1.5, ylim, xlim=NULL,
    factor.names=TRUE, ci.style, 
    type=c("response", "link"), ticks=list(at=NULL, n=5),  
    alternating=TRUE, rotx=0, roty=0, grid=FALSE, layout, rescale.axis=TRUE, 
    transform.x=NULL, ticks.x=NULL,
    key.args=NULL, 
    row=1, col=1, nrow=1, ncol=1, more=FALSE, ...)
{  
    ci.style <- if(missing(ci.style)) NULL else 
        match.arg(ci.style, c("bars", "lines", "none")) 
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
    require(lattice)
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
                        else{ if(ci.style == "lines") {
                            llines(x[good], lower[good], lty=2, col=colors[2])
                            llines(x[good], upper[good], lty=2, col=colors[2])
                        }}
                    }
                    llines(x[good], y[good], lwd=2, col=colors[1], type='b', pch=19, cex=cex, ...)
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
                make.ticks(trans(xlm), link=transform.x[[nm]]$trans, inverse=transform.x[[nm]]$inverse, at=ticks.x$at, n=ticks.x$n)
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
                    llines(x[good], y[good], lwd=2, col=colors[1], ...)
                    if (rug) lrug(x.vals)
                    if (has.se){  
                        if (ci.style == "bars"){
                            larrows(x0=x[good], y0=lower[good], 
                                    x1=x[good], y1=upper[good], 
                                    angle=90, code=3, col=eval(colors[2]), 
                                    length=.125*cex/1.5)
                        }
                        else{ if(ci.style == "lines") {
                            llines(x[good], lower[good], lty=2, col=colors[2])
                            llines(x[good], upper[good], lty=2, col=colors[2])
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
                xlim=trans(xlm),
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
        if (length(which.x) == 0) stop(paste("x.var = '", x.var, "' is not in the model.", sep=""))
        x.var <- which.x
    }
    if (is.character(z.var)) {
        which.z <- which(z.var == predictors)
        if (length(which.z) == 0) stop(paste("z.var = '", z.var, "' is not in the model.", sep=""))
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
            key<-list(title=predictors[z.var], cex.title=1, border=TRUE,
                text=list(as.character(zvals)), 
                lines=list(col=colors[1:length(zvals)], lty=lines[1:length(zvals)], lwd=2), 
                points=list(pch=1:length(zvals)))
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
                        llines(x[sub][good]+os, y[sub][good], lwd=2, type='b', col=colors[i], 
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
                make.ticks(trans(xlm), link=transform.x[[nm]]$trans, inverse=transform.x[[nm]]$inverse, at=ticks.x$at, n=ticks.x$n)
            }
            else {
                trans <- I
                make.ticks(xlm, link=I, inverse=I, at=at, n=n)
            }
            key<-list(title=predictors[z.var], cex.title=1, border=TRUE,
                text=list(as.character(zvals)), 
                lines=list(col=colors[1:length(zvals)], lty=lines[1:length(zvals)], lwd=2))
            key <- c(key, key.args) 
            plot <- xyplot(eval(parse( 
                text=paste("fit ~trans(", predictors[x.var], ")", 
                    if (n.predictors > 2) paste(" |", 
                        paste(predictors[-c(x.var, z.var)])), collapse="*"))),
                strip=function(...) strip.default(..., strip.names=c(factor.names, TRUE)),
# old           panel=function(x, y, subscripts, x.vals, rug, z, ...){
                panel=function(x, y, subscripts, x.vals, rug, z, lower, upper, show.se, ...){
                    if (grid) panel.grid()
                    if (rug) lrug(x.vals)
#new
                    axis.length <- diff(range(x))
#end 
                    for (i in 1:length(zvals)){
                        sub <- z[subscripts] == zvals[i]
                        good <- !is.na(y[sub])
                        llines(x[sub][good], y[sub][good], lwd=2, type='l', col=colors[i], lty=lines[i], cex=cex, ...)
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
                xlim=trans(xlm), 
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
                    else{ if(ci.style == "lines") {
                        llines(x[good], lower[subscripts][good], lty=2, col=colors[2])
                        llines(x[good], upper[subscripts][good], lty=2, col=colors[2])
                    }}
                }
                llines(x[good], y[good], lwd=2, type='b', col=colors[1], pch=19, cex=cex, ...)
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
            make.ticks(trans(xlm), link=transform.x[[nm]]$trans, inverse=transform.x[[nm]]$inverse, at=ticks.x$at, n=ticks.x$n)
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
                llines(x[good], y[good], lwd=2, col=colors[1], ...)
                if (rug) lrug(x.vals)
                if (has.se){  
                    if (ci.style == "bars"){
                            larrows(x0=x[good], y0=lower[subscripts][good], 
                                    x1=x[good], y1=upper[subscripts][good], 
                                    angle=90, code=3, col=eval(colors[2]), 
                                    length=.125*cex/1.5)
                        }
                        else{ if(ci.style == "lines") {
                            llines(x[good], lower[subscripts][good], lty=2, col=colors[2])
                            llines(x[good], upper[subscripts][good], lty=2, col=colors[2])
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
            xlim=trans(xlm),
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


print.effpoly <- function(x, type=c("probability", "logits"), ...){
	type <- match.arg(type)
	x.frame <-as.data.frame(x)
	n.predictors <- length(names(x$x))
	predictors <- names(x.frame)[1:n.predictors]
	y.lev <- x$y.lev
	ylevel.names <- make.names(paste("prob",y.lev))
	colnames(x$prob) <- colnames(x$logit) <- ylevel.names
	y.categories <- matrix(0, nrow=length(x.frame[,predictors[1]]), ncol=length(y.lev))
	for (i in 1:length(y.lev)){
		level <- which(colnames(x$prob)[i] == ylevel.names)
		y.categories[,i] <-  rep(y.lev[level], length(y.categories[,i]))
	}
	y.categories <- as.vector(y.categories)
	y.categories <- factor(y.categories)
	for (i in 1:length(y.lev)){
		cat(paste("\n", gsub(":", "*", x$term), " effect (", type,") for ", y.lev[i], "\n", sep=""))    
		table <- array(if (type == "probability") {x$prob[y.categories==y.lev[i]]}
				else {x$logit[y.categories==y.lev[i]]},     
			dim=sapply(x$variables, function(x) length(x$levels)),
			dimnames=lapply(x$variables, function(x) x$levels))
		print(table, ...)
	}
	if (x$discrepancy > 0.1) cat(paste("\nWarning: There is an average discrepancy of", 
				round(x$discrepancy, 2),
				"percent \n     in the 'safe' predictions for effect", x$term, '\n'))
	invisible(x)
}

summary.effpoly <- function(object, type=c("probability", "logits"), ...){
	type <- match.arg(type)
	x.frame <-as.data.frame(object)
	n.predictors <- length(names(object$x))
	predictors <- names(x.frame)[1:n.predictors]
	y.lev <- object$y.lev
	ylevel.names <- make.names(paste("prob",y.lev))
	colnames(object$prob) <- colnames(object$logit) <- 
		colnames(object$lower.logit) <- colnames(object$upper.logit) <- 
		colnames(object$lower.prob) <- colnames(object$upper.prob)<- ylevel.names
	y.categories <-matrix(0, nrow=length(x.frame[,predictors[1]]), ncol=length(y.lev))
	for (i in 1:length(y.lev)){
		level <- which(colnames(object$prob)[i] == ylevel.names)
		y.categories[,i] <- rep(y.lev[level], length(y.categories[,i]))
	}
	y.categories <- as.vector(y.categories)
	y.categories <- factor(y.categories)
	for (i in 1:length(y.lev)){
		cat(paste("\n", gsub(":", "*", object$term), " effect (" , type, ") for ", y.lev[i], "\n", sep=""))    
		table <- array(if (type == "probability") {object$prob[y.categories==y.lev[i]]}
				else {object$logit[y.categories==y.lev[i]]},     
			dim=sapply(object$variables, function(x) length(x$levels)),
			dimnames=lapply(object$variables, function(x) x$levels))
		print(table, ...)
	}
	if (is.null(object$confidence.level)) return(invisible(NULL))
	for (i in 1:length(y.lev)){
		cat(paste("\n", 'Lower', object$confidence.level*100, 'Percent Confidence Limits for'
				, y.lev[i],'\n'))
		table <- if (type == "probability") object$lower.prob else object$lower.logit
		table <- array(table[y.categories==y.lev[i]],     
			dim=sapply(object$variables, function(x) length(x$levels)),
			dimnames=lapply(object$variables, function(x) x$levels))
		print(table, ...)
	}
	for (i in 1:length(y.lev)){
		cat(paste("\n", 'Upper', object$confidence.level*100, 'Percent Confidence Limits for'
				, y.lev[i],'\n'))
		table <- if (type == "probability") object$upper.prob else object$upper.logit
		table <- array(table[y.categories==y.lev[i]],     
			dim=sapply(object$variables, function(x) length(x$levels)),
			dimnames=lapply(object$variables, function(x) x$levels))
		print(table, ...)
	}
	if (object$discrepancy > 0.1) cat(paste("\nWarning: There is an average discrepancy of", 
				round(object$discrepancy, 2),
				"percent \n     in the 'safe' predictions for effect", object$term, '\n'))
	invisible(NULL)
}

plot.effpoly <- function(x,
    type=c("probability", "logit"),
    x.var=which.max(levels),
    rug=TRUE,
    xlab,
    ylab=paste(x$response, " (", type, ")", sep=""), 
    main=paste(effect, "effect plot"),
    colors, symbols=1:10, lines=1:10, cex=1.5, 
    factor.names=TRUE, ci.style,
    style=c("lines", "stacked"), 
    confint=(style == "lines" && !is.null(x$confidence.level)), 
    transform.x=NULL, ticks.x=NULL, xlim=NULL,
    ylim, rotx=0, alternating=TRUE, roty=0, grid=FALSE,
    layout, key.args=NULL,
    row=1, col=1, nrow=1, ncol=1, more=FALSE, ...){     
    require(lattice)
    ci.style <- if(missing(ci.style)) NULL else 
       match.arg(ci.style, c("bars", "lines", "none"))
    type <- match.arg(type)
    style <- match.arg(style)
    has.se <- !is.null(x$confidence.level) 
    if (confint && !has.se) stop("there are no confidence limits to plot")
    if (style == "stacked"){
        if (type != "probability"){
            type <- "probability"
            warning('type set to "probability" for stacked plot')
        }
        if (confint){
            confint <- FALSE
            warning('confint set to FALSE for stacked plot')
        }
    }
    if (missing(colors)){
        if (style == "stacked"){
            colors <- if (x$model == "multinom") rainbow_hcl(length(x$y.levels))
            else sequential_hcl(length(x$y.levels))
        }
        else colors <- palette()
    }
    effect <- paste(sapply(x$variables, "[[", "name"), collapse="*")
    split <- c(col, row, ncol, nrow)
    n.predictors <- length(names(x$x))
    y.lev <- x$y.lev
    n.y.lev <- length(y.lev)
    ylevel.names <- make.names(paste("prob",y.lev))
    colnames(x$prob) <- colnames(x$logit) <- 
        colnames(x$lower.logit) <- colnames(x$upper.logit) <- 
        colnames(x$lower.prob) <- colnames(x$upper.prob)<- ylevel.names
    x.frame <-as.data.frame(x)
    predictors <- names(x.frame)[1:n.predictors]
    levels <- if (n.predictors==1) length (x.frame[,predictors])
    else sapply(apply(x.frame[,predictors], 2, unique), length)
    if (is.character(x.var)) {
        which.x <- which(x.var == predictors)
        if (length(which.x) == 0) stop(paste("x.var = '", x.var, "' is not in the model.", sep=""))
        x.var <- which.x
    }
    x.vals <- x.frame[, names(x.frame)[x.var]]    
    response <-matrix(0, nrow=nrow(x.frame), ncol=n.y.lev)
    for (i in 1:length(x$y.lev)){
        level <- which(colnames(x$prob)[i] == ylevel.names)
        response[,i] <- rep(x$y.lev[level], length(response[,i]))
    }
    prob <- as.vector(x$prob)
    logit <- as.vector(x$logit)
    response <- as.vector(response)
    if (has.se){
        lower.prob <- as.vector(x$lower.prob)
        upper.prob <- as.vector(x$upper.prob)
        lower.logit <- as.vector(x$lower.logit)
        upper.logit <- as.vector(x$upper.logit)
    }
    response <- factor(response, levels=y.lev)
    data <- data.frame(prob, logit)
    if (has.se) data <- cbind(data, data.frame(lower.prob, upper.prob, lower.logit, upper.logit))
    data[[x$response]] <- response
    for (i in 1:length(predictors)){
        data <-cbind(data, x.frame[predictors[i]])
    }
    levs <- levels(x$data[[predictors[x.var]]])
    n.predictor.cats <- sapply(data[, predictors[-c(x.var)], drop=FALSE], 
        function(x) length(unique(x)))
    if (length(n.predictor.cats) == 0) n.predictor.cats <- 1
    ci.style <- if(is.null(ci.style)) {
       if(is.factor(x$data[[predictors[x.var]]])) "bars" else "lines"} else ci.style
    if( ci.style=="none" ) confint <- FALSE
### no confidence intervals if confint == FALSE or ci.style=="none"
    if (!confint){ # plot without confidence bands
        layout <- if (missing(layout)){
            lay <- c(prod(n.predictor.cats[-(n.predictors - 1)]), 
                prod(n.predictor.cats[(n.predictors - 1)]), 1)
            if (lay[1] > 1) lay else lay[c(2, 1, 3)]
        }
        else layout
        if (style == "lines"){ # line plot
            if (n.y.lev > min(c(length(colors), length(lines), length(symbols))))
                stop(paste('Not enough colors, lines, or symbols to plot', n.y.lev, 'lines'))
            if (is.factor(x$data[[predictors[x.var]]])){ # x-variable a factor
                key <- list(title=x$response, cex.title=1, border=TRUE,
                    text=list(as.character(unique(response))),
                    lines=list(col=colors[1:n.y.lev], lty=lines[1:n.y.lev], lwd=2),
                    points=list(pch=symbols[1:n.y.lev], col=colors[1:n.y.lev]))
                result <- xyplot(eval(if (type=="probability") 
                    parse(text=if (n.predictors==1) 
                        paste("prob ~ as.numeric(", predictors[x.var], ")")
                        else paste("prob ~ as.numeric(", predictors[x.var],") | ", 
                            paste(predictors[-x.var], collapse="*")))
                    else parse(text=if (n.predictors==1) 
                        paste("logit ~ as.numeric(", predictors[x.var], ")")
                        else paste("logit ~ as.numeric(", predictors[x.var],") | ", 
                            paste(predictors[-x.var], collapse="*")))), 
                    strip=function(...) strip.default(..., strip.names=c(factor.names, TRUE)),
                    panel=function(x, y, subscripts, rug, z, x.vals, ...){
                        if (grid) panel.grid()
                        for (i in 1:n.y.lev){
                            sub <- z[subscripts] == y.lev[i]
                            good <- !is.na(y[sub])
                            llines(x[sub][good], y[sub][good], lwd=2, type="b", col=colors[i], lty=lines[i], 
                                pch=symbols[i], cex=cex, ...)
                        }
                    },
                    ylab=ylab,
                    ylim= if (missing(ylim))
                        if (type == "probability") range(prob) else range(logit)
                    else ylim,
                    xlab=if (missing(xlab)) predictors[x.var] else xlab,
                    x.vals=x$data[[predictors[x.var]]], 
                    rug=rug,
                    z=response,
                    scales=list(x=list(at=1:length(levs), labels=levs, rot=rotx),
                        y=list(rot=roty),  
                        alternating=alternating),
                    main=main,
                    key=c(key, key.args),
                    layout=layout,
                    data=data, ...)
                result$split <- split
                result$more <- more
                class(result) <- c("plot.eff", class(result))    		
            }
            else { # x-variable numeric
                nm <- predictors[x.var]
                x.vals <- x$data[[nm]]   
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
                else range.adj(data[nm]) # range(x.vals)
                tickmarks.x <- if ((nm %in% names(transform.x)) && !(is.null(transform.x))){
                    trans <- transform.x[[nm]]$trans
                    make.ticks(trans(xlm), link=transform.x[[nm]]$trans, inverse=transform.x[[nm]]$inverse, at=ticks.x$at, n=ticks.x$n)
                }
                else {
                    trans <- I
                    make.ticks(xlm, link=I, inverse=I, at=at, n=n)
                }
                key <- list(title=x$response, cex.title=1, border=TRUE,
                    text=list(as.character(unique(response))), 
                    lines=list(col=colors[1:n.y.lev], lty=lines[1:n.y.lev], lwd=2))
                result <- xyplot(eval(if (type=="probability") 
                    parse(text=if (n.predictors==1) paste("prob ~ trans(", predictors[x.var], ")")
                        else paste("prob ~ trans(", predictors[x.var],") |", 
                            paste(predictors[-x.var], collapse="*")))
                    else parse(text=if (n.predictors==1) paste("logit ~ trans(", predictors[x.var], ")")
                        else paste("logit ~ trans(", predictors[x.var],") | ", 
                            paste(predictors[-x.var], collapse="*")))), 
                    strip=function(...) strip.default(..., strip.names=c(factor.names, TRUE)),
                    panel=function(x, y, subscripts, rug, z, x.vals, ...){
                        if (grid) panel.grid()
                        if (rug) lrug(x.vals)
                        for (i in 1:n.y.lev){
                            sub <- z[subscripts] == y.lev[i]
                            good <- !is.na(y[sub])
                            llines(x[sub][good], y[sub][good], lwd=2, type="l", col=colors[i], lty=lines[i], ...)
                        }
                    },
                    ylab=ylab,
                    xlim=trans(xlm),
                    ylim= if (missing(ylim))
                        if (type == "probability") range(prob) else range(logit)
                    else ylim,
                    xlab=if (missing(xlab)) predictors[x.var] else xlab,
                    x.vals=x$data[[predictors[x.var]]], 
                    rug=rug,
                    z=response,
                    scales=list(x=list(at=tickmarks.x$at, labels=tickmarks.x$labels, rot=rotx), y=list(rot=roty),
                        alternating=alternating),
                    main=main,
                    key=c(key, key.args),
                    layout=layout,
                    data=data, ...)
                result$split <- split
                result$more <- more
                class(result) <- c("plot.eff", class(result))			
            }
        }
        else { # stacked plot
            if (n.y.lev > length(colors))
                stop(paste('Not enough colors to plot', n.y.lev, 'regions'))
            key <- list(text=list(lab=rev(y.lev)), rectangle=list(col=rev(colors[1:n.y.lev])))
            if (is.factor(x$data[[predictors[x.var]]])){ # x-variable a factor
                result <- barchart(eval(parse(text=if (n.predictors == 1) 
                    paste("prob ~ ", predictors[x.var], sep="")
                    else paste("prob ~ ", predictors[x.var]," | ", 
                        paste(predictors[-x.var], collapse="*")))), 
                    strip=function(...) strip.default(..., strip.names=c(factor.names, TRUE)),
                    groups = response,
                    col=colors,
                    horizontal=FALSE, 
                    stack=TRUE, 
                    data=data, 
                    ylim=if (missing(ylim)) 0:1 else ylim,
                    ylab=ylab, 
                    xlab=if (missing(xlab)) predictors[x.var] else xlab,
                    scales=list(x=list(rot=rotx), y=list(rot=roty), 
                        alternating=alternating),
                    main=main,
                    key=c(key, key.args),
                    layout=layout)
                result$split <- split
                result$more <- more
                class(result) <- c("plot.eff", class(result))			
            }
            else { # x-variable numeric
                nm <- predictors[x.var]
                x.vals <- x$data[[nm]]   
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
                else range.adj(data[nm]) # range(x.vals)
                tickmarks.x <- if ((nm %in% names(transform.x)) && !(is.null(transform.x))){
                    trans <- transform.x[[nm]]$trans
                    make.ticks(trans(xlm), link=transform.x[[nm]]$trans, inverse=transform.x[[nm]]$inverse, at=ticks.x$at, n=ticks.x$n)
                }
                else {
                    trans <- I
                    make.ticks(xlm, link=I, inverse=I, at=at, n=n)
                }
                result <- densityplot(eval(parse(text=if (n.predictors == 1)
                    paste("~ trans(", predictors[x.var], ")", sep="")
                    else paste("~ trans(", predictors[x.var], ") | ",
                        paste(predictors[-x.var], collapse="*")))),
                    probs=x$prob,
                    strip=function(...) strip.default(..., strip.names=c(factor.names, TRUE)),
                    panel =  function(x, subscripts, rug, x.vals, probs=probs, col=colors, ...){
                        fill <- function(x, y1, y2, col){
                            if (length(y2) == 1) y2 <- rep(y2, length(y1))
                            if (length(y1) == 1) y1 <- rep(y1, length(y2))
                            panel.polygon(c(x, rev(x)), c(y1, rev(y2)), col=col)
                        }
                        n <- ncol(probs)
                        Y <- t(apply(probs[subscripts,], 1, cumsum))
                        fill(x, 0, Y[,1], col=col[1])
                        for (i in 2:n){
                            fill(x, Y[,i-1], Y[,i], col=col[i])
                        }
                        if (rug) lrug(x.vals)
                    },
                    rug=rug,
                    x.vals=x$data[[predictors[x.var]]],
                    data=x$x,
                    xlim=trans(xlm),
                    ylim=if (missing(ylim)) 0:1 else ylim,
                    ylab=ylab,
                    xlab=if (missing(xlab)) predictors[x.var] else xlab,
                    scales=list(x=list(at=tickmarks.x$at, labels=tickmarks.x$labels, rot=rotx), y=list(rot=roty),
                        alternating=alternating),
                    main=main,
                    key=c(key, key.args),
                    layout=layout, ...)
                result$split <- split
                result$more <- more
                class(result) <- c("plot.eff", class(result))
            }
        }
    }
### with confidence banks
    else{ # plot with confidence bands
        layout <- if(missing(layout)) c(prod(n.predictor.cats), length(levels(response)), 1) 
        else layout
        if (type == "probability"){
            lower <- lower.prob
            upper <- upper.prob
        }
        else {
            lower <- lower.logit
            upper <- upper.logit
        }
### factor
        if (is.factor(x$data[[predictors[x.var]]])){ # x-variable a factor
            levs <- levels(x$data[[predictors[x.var]]])
            result <- xyplot(eval(if (type=="probability") 
                parse(text=if (n.predictors==1) 
                    paste("prob ~ as.numeric(", predictors[x.var],") |", x$response)
                    else paste("prob ~ as.numeric(", predictors[x.var],") |", 
                        paste(predictors[-x.var], collapse="*"), 
                        paste("*", x$response)))
                else parse(text=if (n.predictors==1) 
                    paste("logit ~ as.numeric(", predictors[x.var],") |", x$response)
                    else paste("logit ~ as.numeric(", predictors[x.var],")|", 
                        paste(predictors[-x.var], collapse="*"), 
                        paste("*", x$response)))),
                par.strip.text=list(cex=0.8),							
                strip=function(...) strip.default(..., strip.names=c(factor.names, TRUE)),
                panel=function(x, y, subscripts, x.vals, rug, lower, upper, ... ){
                    if (grid) panel.grid()
                    good <- !is.na(y)
                    llines(x[good], y[good], lwd=2, type="b", pch=19, col=colors[1], cex=cex, ...)
                    if (ci.style == "bars"){
                        larrows(x0=x[good], y0=lower[subscripts+as.numeric(rownames(data)[1])-1][good], 
                            x1=x[good], y1=upper[subscripts+as.numeric(rownames(data)[1])-1][good], 
                            angle=90, code=3, col=colors[2], length=0.125*cex/1.5)
                    }
                    else { if(ci.style == "lines"){
                        llines(x[good], lower[subscripts+as.numeric(rownames(data)[1])-1][good], lty=2, col=colors[2])
                        llines(x[good], upper[subscripts+as.numeric(rownames(data)[1])-1][good], lty=2, col=colors[2])
                    } }
                },
                ylab=ylab,
                ylim= if (missing(ylim)) c(min(lower), max(upper)) else ylim,
                xlab=if (missing(xlab)) predictors[x.var] else xlab,
                main=main,
                x.vals=x$data[[predictors[x.var]]],
                rug=rug,
                lower=lower,
                upper=upper, 
                scales=list(x=list(at=1:length(levs), labels=levs, rot=rotx), 
                    y=list(rot=roty), alternating=alternating),
                layout=layout,
                data=data, ...)
            result$split <- split
            result$more <- more
            class(result) <- c("plot.eff", class(result))
        }
        else { # x-variable numeric
            nm <- predictors[x.var]
            x.vals <- x$data[[nm]]   
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
            else range.adj(data[nm]) # range(x.vals)
            tickmarks.x <- if ((nm %in% names(transform.x)) && !(is.null(transform.x))){
                trans <- transform.x[[nm]]$trans
                make.ticks(trans(xlm), link=transform.x[[nm]]$trans, inverse=transform.x[[nm]]$inverse, at=ticks.x$at, n=ticks.x$n)
            }
            else {
                trans <- I
                make.ticks(xlm, link=I, inverse=I, at=at, n=n)
            }
            result <- xyplot(eval(if (type=="probability") 
                parse(text=if (n.predictors==1) 
                    paste("prob ~ trans(", predictors[x.var],") |", x$response)
                    else paste("prob ~ trans(", predictors[x.var],") |", 
                        paste(predictors[-x.var], collapse="*"), 
                        paste("*", x$response)))
                else parse(text=if (n.predictors==1) 
                    paste("logit ~ trans(", predictors[x.var],") |", x$response)
                    else paste("logit ~ trans(", predictors[x.var],") |", 
                        paste(predictors[-x.var], collapse="*"), 
                        paste("*", x$response)))
            ),
                par.strip.text=list(cex=0.8),							
                strip=function(...) strip.default(..., strip.names=c(factor.names, TRUE)),
                panel=function(x, y, subscripts, x.vals, rug, lower, upper, ... ){
                    if (grid) panel.grid()
                    if (rug) lrug(x.vals)
                    good <- !is.na(y)
                    llines(x[good], y[good], lwd=2, col=colors[1], ...)
                    if (ci.style == "bars"){
                        larrows(x0=x[good], y0=lower[subscripts+as.numeric(rownames(data)[1])-1][good], 
                            x1=x[good], y1=upper[subscripts+as.numeric(rownames(data)[1])-1][good], 
                            angle=90, code=3, col=colors[2], length=0.125*cex/1.5)
                    }
                    else { if(ci.style == "lines"){
                        llines(x[good], lower[subscripts+as.numeric(rownames(data)[1])-1][good], lty=2, col=colors[2])
                        llines(x[good], upper[subscripts+as.numeric(rownames(data)[1])-1][good], lty=2, col=colors[2])
                    } }
                },
                ylab=ylab,
                xlim=trans(xlm),
                ylim= if (missing(ylim)) c(min(lower), max(upper)) else ylim,
                xlab=if (missing(xlab)) predictors[x.var] else xlab,
                main=main,
                x.vals=x$data[[predictors[x.var]]],
                rug=rug,
                lower=lower,
                upper=upper, 
                scales=list(y=list(rot=roty), x=list(at=tickmarks.x$at, labels=tickmarks.x$labels, rot=rotx),
                    alternating=alternating),
                layout=layout,
                data=data, ...)
            result$split <- split
            result$more <- more
            class(result) <- c("plot.eff", class(result))
        }
    }
    result
}

print.efflatent <- function(x, ...){
	cat(paste("\n", gsub(":", "*", x$term), 'effect\n'))
	table <- array(x$fit,     
		dim=sapply(x$variables, function(x) length(x$levels)),
		dimnames=lapply(x$variables, function(x) x$levels))
	print(table, ...)
	cat("\nThresholds:\n")
	print(x$thresholds, ...)
	if (x$discrepancy > 0.1) cat(paste("\nWarning: There is an average discrepancy of", 
				round(x$discrepancy, 3),
				"percent \n     in the 'safe' predictions for effect", x$term, '\n'))
	invisible(x)
}

summary.efflatent <- function(object, ...){
	result <- list()
	result$header <- paste("\n", gsub(":", "*", object$term), 'effect\n')
	result$effect <- array(object$fit,     
		dim=sapply(object$variables, function(x) length(x$levels)),
		dimnames=lapply(object$variables, function(x) x$levels))
	if (!is.null(object$se)){
		result$lower.header <- paste('\n Lower', round(100*object$confidence.level, 2), 
			'Percent Confidence Limits\n')
		result$lower <- array(object$lower,   
			dim=sapply(object$variables, function(x) length(x$levels)),
			dimnames=lapply(object$variables, function(x) x$levels))
		result$upper.header <- paste('\n Upper', round(100*object$confidence.level, 2),
			'Percent Confidence Limits\n')
		result$upper <- array(object$upper,   
			dim=sapply(object$variables, function(x) length(x$levels)),
			dimnames=lapply(object$variables, function(x) x$levels))
	}
	result$thresholds <- object$thresholds
	if (object$discrepancy > 0.1) result$warning <- paste("\nWarning: There is an average discrepancy of", 
			round(object$discrepancy, 3),
			"percent \n     in the 'safe' predictions for effect", object$term, '\n')
	class(result) <- "summary.eff"
	result
}
