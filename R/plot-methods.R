
# plot.eff method for effects package, moved here from plot-summary-print-methods.R
# The plot.effpoly method remains there for now.
# 2013-10-17: Added use.splines keyword to plot.eff. Sandy
# 2013-10-17: Made ci.style="bands" default for variates; allow "bands" if multiline=TRUE
# 2013-10-29: fixed plot.eff() to handle factors with "valid" NA level. J. Fox
# 2014-03-03: modified plot.eff() to handle partial residuals. J. Fox
# 2014-09-20: fixed plot.eff() to work with partial residuals when rescale.axis=FALSE;
#             added smooth.residuals argument. J. Fox
# 2014-10-10: namespace fixes. J. Fox
# 2014-12-05: made key.args more flexible. J. Fox
# 2015-03-22: use wide columns by default only when x for legend not set. J. Fox
# 2015-03-25: use non-robust loess smooth for partial residuals for non-Gaussian families. J. Fox
# 2015-03-25: rationalized type and rescale.axis args to plot.eff(); deprecated rescale.axis arg. J. Fox
# 2015-05-28: added residuals.smooth.color argument. J. Fox
# 2015-08-28: added residuals.cex argument. J. Fox
# 2016-03-01: move computation of partial residuals to the plot.eff() method. J. Fox
# 2016-05-22: modified make.ticks() to avoid possible failure due to floating-point inaccuracy. J. Fox
# 2016-08-31: fixed plotting with partial residuals with various scalings of y-axis and x-axis. J. Fox
# 2016-09-16: added show.strip.values argument to plot.eff(). J. Fox
# 2017-06-12: fixed bug in plot.eff() for multiline displays with many conditioning variables. J. Fox
# 2017-07-15: modified plot.eff() to consolidate arguments and use lattice theme. J. Fox
# 2017-08-09: small bug fixes, reorganized axes=list(x=list()) argument. J. Fox
# 2017-08-17: tweaked layout. J. Fox
# 2017-08-23: Fixed bug with the lattice=list(array()) argument in plot.efflist --- lattice was as
#              an argument to the next method twice
# 2017-08-23: plot.eff, in key.args, set default for between.columns=0
# 2017-08-20: reintroduce legacy arguments for plot.eff()
# 2017-09-10: use replacement for grid.panel()

# the following functions aren't exported

find.legend.columns <- function(n, target=min(4, n)){
  rem <- n %% target
  if (rem != 0 && rem < target/2) target <- target - 1
  target
}


make.ticks <- function(range, link, inverse, at, n) {
  warn <- options(warn=-1)
  on.exit(options(warn))
  link <- if (is.null(link))
    function(x) nlm(function(y) (inverse(y) - x)^2,
                    mean(range))$estimate
  else link
  if (is.null(n)) n <- 5
  labels <- if (is.null(at)){
    range.labels <- sapply(range, inverse)
    labels <- grid.pretty(range.labels)
  }
  else at
  ticks <- try(sapply(labels, link), silent=TRUE)
  if (inherits(ticks, "try-error")){
    ticks <- seq(range[1], range[2], length=n)
  }
  list(at=ticks, labels=format(labels))
}

range.adj <- function(x){
  range <- range(x, na.rm=TRUE)
  c(range[1] - .025*(range[2] - range[1]),
    range[2] + .025*(range[2] - range[1]))
}

# added, modified from http://www.r-bloggers.com/confidence-bands-with-lattice-and-r/

panel.bands <- function(x, y, upper, lower, fill, col,
                        subscripts, ..., font, fontface, use.splines=FALSE)
{
  if(!missing(subscripts)) {
    upper <- upper[subscripts]
    lower <- lower[subscripts]
  }
  if (use.splines){
    up <- spline(x, upper)
    down <- spline(x, lower)
    x <- up$x
    upper <- up$y
    lower <- down$y
  }
  panel.polygon(c(x, rev(x)), c(upper, rev(lower)),
                col = fill, fill=fill, border = FALSE,
                ...)
}


# modified by Michael Friendly: added key.args:
# modified by Michael Friendly: added ci.style="bands"
# modified by Michael Friendly: added lwd= argument for llines (not used elsewhere)
# modified by Michael Friendly: added alpha.band= argument for ci.style="bands"

spline.llines <- function(x, y, ...) llines(spline(x, y), ...)

plot.eff <- function(x, x.var, z.var=which.min(levels), main=paste(effect, "effect plot"),
                     symbols=TRUE, lines=TRUE, axes, confint, partial.residuals, id, lattice,
                     ...,
                     # legacy arguments:
                     multiline, rug, xlab, ylab, colors, cex, lty, lwd, ylim, xlim, factor.names, ci.style,
                     band.transparency, band.colors, type, ticks, alternating, rotx, roty, grid, layout,
                     rescale.axis, transform.x, ticks.x, show.strip.values, key.args, use.splines,
                     residuals.color, residuals.pch, residuals.cex, smooth.residuals,
                     residuals.smooth.color, show.fitted, span)
{
  closest <- function(x, x0) apply(outer(x, x0, FUN=function(x, x0) abs(x - x0)), 1, which.min)
  .mod <- function(a, b) ifelse( (d <- a %% b) == 0, b, d)
  .modc <- function(a) .mod(a, length(colors))
  .mods <- function(a) .mod(a, length(symbols))
  .modl <- function(a) .mod(a, length(lines))
  .modb <- function(a) .mod(a, length(band.colors))
  
  if (!is.logical(lines) && !is.list(lines)) lines <- list(lty=lines)
  lines <- applyDefaults(lines,
                         defaults=list(multiline=is.null(x$se), lty=trellis.par.get("superpose.line")$lty,
                                       lwd=trellis.par.get("superpose.line")$lwd[1], col=trellis.par.get("superpose.line")$col, splines=TRUE),
                         onFALSE=list(multiline=FALSE, lty=0, lwd=0, col=rgb(1, 1, 1, alpha=0), splines=FALSE),
                         arg="lines")
  if (missing(multiline)) multiline <- lines$multiline
  if (missing(lwd)) lwd <- lines$lwd
  if (missing(colors)) colors <- lines$col
  if (missing(use.splines)) use.splines <- lines$splines
  lines <- if (missing(lty)) lines$lty else lty
  
  if (!is.logical(symbols) && !is.list(symbols)) symbols <- list(pch=symbols)
  symbols <- applyDefaults(symbols,
                           defaults=list(pch=trellis.par.get("superpose.symbol")$pch, cex=trellis.par.get("superpose.symbol")$cex[1]),
                           onFALSE=list(pch=NA_integer_, cex=0),
                           arg="symbols")
  cex <- symbols$cex
  symbols <- symbols$pch
  
  if (missing(axes)) axes <- NULL
  axes <- applyDefaults(axes, defaults=list(
    x=list(rotate=0, rug=TRUE),
    y=list(lab=NA, lim=NA, ticks=list(at=NULL, n=5), type="rescale", rotate=0),
    alternating=TRUE, grid=FALSE),
    arg="axes")
  x.args <- applyDefaults(axes$x, defaults=list(rotate=0, rug=TRUE), arg="axes$x")
  
  if (missing(xlab)) {
    xlab.arg <- FALSE
    xlab <- list()
  }
  if (missing(xlim)) {
    xlim.arg <- FALSE
    xlim <- list()
  }
  if (missing(ticks.x)) {
    ticks.x.arg <- FALSE
    ticks.x <- list()
  }
  if (missing(transform.x)) {
    transform.x.arg <- FALSE
    transform.x <- list()
  }
  if (missing(rotx)) rotx <- x.args$rotate
  if (missing(rug)) rug <- x.args$rug
  x.args$rotate <- NULL
  x.args$rug <- NULL
  x.pred.names <- names(x.args)
  if (length(x.pred.names) > 0){
    for (pred.name in x.pred.names){
      x.pred.args <- applyDefaults(x.args[[pred.name]],
                                   defaults=list(lab=NULL, lim=NULL, ticks=NULL, transform=NULL),
                                   arg=paste0("axes$x$", pred.name))
      if (!xlab.arg) xlab[[pred.name]] <- x.pred.args$lab
      if (!xlim.arg) xlim[[pred.name]] <- x.pred.args$lim
      if (!ticks.x.arg) ticks.x[[pred.name]] <- x.pred.args$ticks
      if (!transform.x.arg) transform.x[[pred.name]] <- x.pred.args$transform
    }
  }
  if (length(xlab) == 0) xlab <- NA
  if (length(xlim) == 0) xlim <- NA
  if (length(ticks.x) == 0) ticks.x <- NA
  if (length(transform.x) == 0) transform.x <- NA
  
  y.args <- applyDefaults(axes$y, defaults=list(lab=NA, lim=NA, ticks=list(at=NULL, n=5), type="rescale", rotate=0), arg="axes$y")
  if (missing(ylab)) ylab <- y.args$lab
  if (missing(ylim)) ylim <- y.args$lim
  if (missing(ticks)) ticks <- y.args$ticks
  if (missing(type)) type <- y.args$type
  if (!missing(rescale.axis)) type <- if (rescale.axis) "rescale" else "response"
  type <- match.arg(type, c("rescale", "response", "link"))
  if (missing(roty)) roty <- y.args$rotate
  if (missing(alternating)) alternating <- axes$alternating
  if (missing(grid)) grid <- axes$grid
  
  if (missing(confint) || isTRUE(confint)) confint <- NULL
  confint <- applyDefaults(confint,
                           defaults=list(style=NULL, alpha=0.15, col=colors),
                           onFALSE=list(style="none", alpha=0, col=NA_integer_),
                           arg="confint")
  if (missing(ci.style)) ci.style <- confint$style
  if (missing(band.transparency)) band.transparency <- confint$alpha
  if (missing(band.colors)) band.colors <- confint$col
  if(!is.null(ci.style)) ci.style <- match.arg(ci.style, c("auto", "bars", "lines", "bands", "none"))
  
  if (missing(partial.residuals)) partial.residuals <- NULL
  if (is.logical(partial.residuals)) partial.residuals <- list(plot=partial.residuals)
  partial.residuals <- applyDefaults(partial.residuals, defaults=list(
    plot=!is.null(x$residuals), fitted=FALSE, col=colors[2], pch=1, cex=1, smooth=TRUE, 
    span=2/3, smooth.col=colors[2], lty=lines[1], lwd=lwd),
    arg="partial.residuals")
  if (missing(show.fitted)) show.fitted <- partial.residuals$fitted
  if (missing(residuals.color)) residuals.color <- partial.residuals$col
  if (missing(residuals.pch)) residuals.pch <- partial.residuals$pch
  if (missing(residuals.cex)) residuals.cex <- partial.residuals$cex
  if (missing(smooth.residuals)) smooth.residuals <- partial.residuals$smooth
  if (missing(residuals.smooth.color)) residuals.smooth.color <- partial.residuals$smooth.col
  residuals.lty <- partial.residuals$lty
  residuals.lwd <- partial.residuals$lwd
  if (missing(span)) span <- partial.residuals$span
  partial.residuals <- partial.residuals$plot
  
  if (missing(id) || isFALSE(id)) {
    id.n <- 0
    id.cex <- 0
    id.col <- NULL
    id.labels <- NULL
  }
  else {
    id <- applyDefaults(id, list(
      n=2, cex=0.75, col=residuals.color, labels=NULL
    ), arg="id")
    id.n <- id$n
    id.col <- id$col
    id.cex <- id$cex
    id.labels <- id$labels
  }
  
  if (missing(lattice)) lattice <- NULL
  lattice <- applyDefaults(lattice, defaults=list(
    layout=NULL, key.args=list(),
    strip=list(factor.names=TRUE, values=!partial.residuals),
    array=list(row=1, col=1, nrow=1, ncol=1, more=FALSE),
    arg="lattice"
  ))
  if (missing(layout)) layout <- lattice$layout
  if (missing(key.args)){
    lattice$key.args[["between.columns"]] <- if(is.null(lattice$key.args[["between.columns"]])) 0 else
      lattice$key.args[["between.columns"]]
    key.args <- lattice$key.args
  }
  strip.args <- applyDefaults(lattice$strip, defaults=list(factor.names=TRUE, values=!partial.residuals), arg="lattice$strip")
  if (missing(factor.names)) factor.names <- strip.args$factor.names
  if (missing(show.strip.values)) show.strip.values <- strip.args$values
  array.args <- applyDefaults(lattice$array, defaults=list(row=1, col=1, nrow=1, ncol=1, more=FALSE), arg="lattice$array")
  row <- array.args$row
  col <- array.args$col
  nrow <- array.args$nrow
  ncol <- array.args$ncol
  more <- array.args$more
  
  if (smooth.residuals && !is.null(x$family)){
    loess.family <- if (x$family == "gaussian") "symmetric" else "gaussian"
  }
  
  switch(type,
         rescale = {
           type <- "response"
           rescale.axis <- TRUE
         },
         response = {
           type <- "response"
           rescale.axis <- FALSE
         },
         link = {
           type <- "link"
           rescale.axis <- TRUE
         }
  )
  levels <- sapply(x$variables, function(z) length(as.vector(z[["levels"]])))
  thresholds <- x$thresholds
  has.thresholds <- !is.null(thresholds)
  effect.llines <- llines
  if (is.na(ylab)){
    ylab <- if (has.thresholds) paste(x$response, ": ", paste(x$y.levels, collapse=", "), sep="")
    else x$response
  }
  if (has.thresholds){
    threshold.labels <- abbreviate(x$y.levels, minlength=1)
    threshold.labels <- paste(" ",
                              paste(threshold.labels[-length(threshold.labels)], threshold.labels[-1], sep=" - "),
                              " ", sep="")
  }
  original.link <- trans.link <- x$transformation$link
  original.inverse <- trans.inverse <- x$transformation$inverse
  residuals <- if (partial.residuals) x$residuals else NULL
  if (!is.null(residuals) && !is.null(id.labels)) names(residuals) <- id.labels
  partial.residuals.range <- x$partial.residuals.range
  
  if (!rescale.axis){
    x$lower[!is.na(x$lower)] <- trans.inverse(x$lower[!is.na(x$lower)])
    x$upper[!is.na(x$upper)] <- trans.inverse(x$upper[!is.na(x$upper)])
    x$fit[!is.na(x$fit)] <- trans.inverse(x$fit)[!is.na(x$fit)]
    trans.link <- trans.inverse <- I
  }
  x.all <- x$x.all
  split <- c(col, row, ncol, nrow)
  if (missing(x.var)) x.var <- x$x.var
  if (!is.null(x.var) && is.numeric(x.var)) x.var <- names(x.var)
  x.data <- x$data
  effect <- paste(sapply(x$variables, "[[", "name"), collapse="*")
  vars <- x$variables
  x <- as.data.frame(x, transform=I)
  for (i in 1:length(vars)){
    if (!(vars[[i]]$is.factor)) next
    x[,i] <- factor(x[,i], levels=vars[[i]]$levels, exclude=NULL)
  }
  has.se <- !is.null(x$se)
  n.predictors <- ncol(x) - 1 - 3*has.se
  if (n.predictors == 1){
    predictor <- names(x)[1]
    if (is.list(xlab)) xlab <- xlab[[predictor]]
    ### factor no other predictors
    if (is.factor(x[,1])){
      ci.style <- if(is.null(ci.style) || ci.style == "auto") "bars" else ci.style
      range <- if(has.se & ci.style!="none")
        range(c(x$lower, x$upper), na.rm=TRUE) else range(x$fit, na.rm=TRUE)
      ylim <- if (!any(is.na(ylim))) ylim else c(range[1] - .025*(range[2] - range[1]),
                                                 range[2] + .025*(range[2] - range[1]))
      tickmarks <- if (type == "response" && rescale.axis) make.ticks(ylim,
                                                                      link=trans.link, inverse=trans.inverse, at=ticks$at, n=ticks$n)
      else make.ticks(ylim, link=I, inverse=I, at=ticks$at, n=ticks$n)
      levs <- levels(x[,1])
      plot <- xyplot(eval(parse(
        text=paste("fit ~ as.numeric(", names(x)[1], ")"))),
        strip=function(...) strip.default(..., strip.names=c(factor.names, TRUE)),
        panel=function(x, y, lower, upper, has.se, ...){
          if (grid) ticksGrid(x=1:length(levs), y=tickmarks$at)
          good <- !is.na(y)
          if (has.se){
            if (ci.style == "bars"){
              larrows(x0=x[good], y0=lower[good], x1=x[good], y1=upper[good], angle=90,
                      code=3, col=colors[.modc(2)], length=0.125*cex/1.5)
            }
            else if(ci.style == "lines") {
              effect.llines(x[good], lower[good], lty=2, col=colors[.modc(2)])
              effect.llines(x[good], upper[good], lty=2, col=colors[.modc(2)])
            }
            else{ if(ci.style == "bands") {
              panel.bands(x[good], y[good], upper[good], lower[good], fill=band.colors[1],
                          alpha=band.transparency, use.splines=FALSE)
            }}
          }
          effect.llines(x[good], y[good], lwd=lwd, col=colors[1], lty=lines, type='b', pch=19, cex=cex, ...)
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
        xlab=if (is.na(xlab)) names(x)[1] else xlab,
        scales=list(x=list(at=1:length(levs), labels=levs, rot=rotx),
                    y=list(at=tickmarks$at, labels=tickmarks$labels, rot=roty),
                    alternating=alternating, y=roty),
        main=main,
        lower=x$lower, upper=x$upper, has.se=has.se, data=x, ...)
      result <- update(plot, layout = if (is.null(layout)) c(0, prod(dim(plot)))
                       else layout)
      result$split <- split
      result$more <- more
      class(result) <- c("plot.eff", class(result))
    }
    ### variate, no other predictors  ***
    else {
      effect.llines <- if(use.splines) spline.llines else effect.llines
      ci.style <- if(is.null(ci.style) || ci.style == "auto") "bands" else ci.style
      range <- if(has.se && ci.style!="none")
        range(c(x$lower, x$upper), na.rm=TRUE) else range(x$fit, na.rm=TRUE)
      
      ylim <- if (!any(is.na(ylim))) ylim
      else if (is.null(residuals)) c(range[1] - .025*(range[2] - range[1]), range[2] + .025*(range[2] - range[1]))
      else if (rescale.axis) c(min(partial.residuals.range[1], range[1] - .025*(range[2] - range[1])),
                               max(partial.residuals.range[2], range[2] + .025*(range[2] - range[1])))
      else c(min(original.inverse(partial.residuals.range[1]), range[1] - .025*(range[2] - range[1])),
             max(original.inverse(partial.residuals.range[2]), range[2] + .025*(range[2] - range[1])))
      tickmarks <- if (type == "response" && rescale.axis) make.ticks(ylim,
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
        make.ticks(trans(xlm), link=transform.x[[nm]]$trans, inverse=transform.x[[nm]]$inverse, at=at, n=n)
      }
      else {
        trans <- I
        make.ticks(xlm, link=I, inverse=I, at=at, n=n)
      }
      
      if (is.null(x.var)){
        if (!is.null(residuals)){
          x.var <- names(x)[1]
        }
        else x.var <-  which.max(levels)
      }
      if (!is.null(residuals)) x.fit <- x.data[, predictor]
      if (is.numeric(x.var)) x.var <- predictor
      plot <- xyplot(eval(parse(
        text=paste("fit ~ trans(", x.var, ")"))),
        strip=function(...) strip.default(..., strip.names=c(factor.names, TRUE)),
        panel=function(x, y, x.vals, rug, lower, upper, has.se, ...){
          if (grid) ticksGrid(x=tickmarks.x$at, y=tickmarks$at)
          good <- !is.na(y)
          axis.length <- diff(range(x))
          effect.llines(x[good], y[good], lwd=lwd, col=colors[1], ...)
          if (rug && is.null(residuals)) lrug(trans(x.vals))
          if (has.se){
            if (ci.style == "bars"){
              larrows(x0=x[good], y0=lower[good],
                      x1=x[good], y1=upper[good],
                      angle=90, code=3, col=eval(colors[.modc(2)]),
                      length=.125*cex/1.5)
            }
            else if(ci.style == "lines") {
              effect.llines(x[good], lower[good], lty=2, col=colors[.modc(2)])
              effect.llines(x[good], upper[good], lty=2, col=colors[.modc(2)])
            }
            else{ if(ci.style == "bands") {
              panel.bands(x[good], y[good], upper[good], lower[good], fill=band.colors[1],
                          alpha=band.transparency, use.splines=use.splines)
            }}
          }
          if (has.thresholds){
            panel.abline(h=thresholds, lty=3)
            panel.text(rep(current.panel.limits()$xlim[1], length(thresholds)),
                       thresholds, threshold.labels, adj=c(0,0), cex=0.75)
            panel.text(rep(current.panel.limits()$xlim[2], length(thresholds)),
                       thresholds, threshold.labels, adj=c(1,0), cex=0.75)
          }
          if (!is.null(residuals)){
            fitted <- y[good][closest(trans(x.fit), x[good])]
            partial.res <- if (!rescale.axis) original.inverse(original.link(fitted) + residuals)
            else fitted + residuals
            lpoints(trans(x.fit), partial.res, col=residuals.color, pch=residuals.pch, cex=residuals.cex)
            if (show.fitted) lpoints(trans(x.fit), fitted, pch=16, col=residuals.color)  # REMOVE ME
            if (smooth.residuals){
              llines(loess.smooth(trans(x.fit), partial.res, span=span, family=loess.family), lwd=residuals.lwd, lty=residuals.lty, col=residuals.smooth.color)
            }
            if (id.n > 0){
              M <- cbind(trans(x.fit), partial.res)
              md <- mahalanobis(M, colMeans(M), cov(M))
              biggest <- order(md, decreasing=TRUE)[1:id.n]
              pos <- ifelse(trans(x.fit[biggest]) > mean(current.panel.limits()$xlim), 2, 4)
              ltext(trans(x.fit[biggest]), partial.res[biggest], 
                    names(partial.res)[biggest], pos=pos, col=id.col, cex=id.cex)
            }
          }
          
        },
        ylim=ylim,
        xlim=suppressWarnings(trans(xlm)),
        ylab=ylab,
        xlab=if (is.na(xlab)) names(x)[1] else xlab,
        x.vals=x.vals, rug=rug,
        main=main,
        lower=x$lower, upper=x$upper, has.se=has.se, data=x,
        scales=list(y=list(at=tickmarks$at, labels=tickmarks$labels, rot=roty),
                    x=list(at=tickmarks.x$at, labels=tickmarks.x$labels, rot=rotx), alternating=alternating), ...)
      result <- update(plot, layout = if (is.null(layout)) c(0, prod(dim(plot)))
                       else layout)
      result$split <- split
      result$more <- more
      class(result) <- c("plot.eff", class(result))
    }
    return(result)
  }
  ###  more than one predictor
  predictors <- names(x)[1:n.predictors]
  levels <- sapply(apply(x[,predictors], 2, unique), length)
  if (is.null(x.var)){
    if (!is.null(residuals)){
      x.var <- names(x)[1]
    }
    else x.var <-  which.max(levels)
  }
  if (is.list(xlab)) xlab <- xlab[[x.var]]
  if (!is.null(residuals)) x.fit <- x.data[, x.var]
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
    if (!is.null(residuals)) warning("partial residuals are not displayed in a multiline plot")
    ci.style <- if(is.null(ci.style)) "none" else ci.style
    if(ci.style == "lines") {
      cat("Confidence interval style 'lines' changed to 'bars'\n")
      ci.style <- "bars"}
    range <- if (has.se && ci.style !="none")
      range(c(x$lower, x$upper), na.rm=TRUE) else range(x$fit, na.rm=TRUE)
    ylim <- if (!any(is.na(ylim))) ylim else c(range[1] - .025*(range[2] - range[1]),
                                               range[2] + .025*(range[2] - range[1]))
    tickmarks <- if (type == "response" && rescale.axis) make.ticks(ylim, link=trans.link,
                                                                    inverse=trans.inverse, at=ticks$at, n=ticks$n)
    else make.ticks(ylim, link=I, inverse=I, at=ticks$at, n=ticks$n)
    zvals <- unique(x[, z.var])
    ### multiline factor
    if (is.factor(x[,x.var])){
      if (ci.style == "auto") ci.style <- "bars"
      levs <- levels(x[,x.var])
      key <- list(title=predictors[z.var], cex.title=1, border=TRUE,
                  text=list(as.character(zvals)),
                  lines=list(col=colors[.modc(1:length(zvals))], lty=lines[.modl(1:length(zvals))], lwd=lwd),
                  points=list(col=colors[.modc(1:length(zvals))], pch=symbols[.mods(1:length(zvals))]),
                  columns = if ("x" %in% names(key.args)) 1 else find.legend.columns(length(zvals)))
      for (k in names(key.args)) key[k] <- key.args[k]
      if (show.strip.values && n.predictors > 2){
        for (pred in predictors[-c(x.var, z.var)]){
          x[[pred]] <- as.factor(x[[pred]])
        }
      }
      plot <- xyplot(eval(parse(
        text=paste("fit ~ as.numeric(", predictors[x.var], ")",
                   if (n.predictors > 2) paste(" |",
                                               paste(predictors[-c(x.var, z.var)], collapse="*"))))),
        strip=function(...) strip.default(..., strip.names=c(factor.names, TRUE), sep=" = "),
        panel=function(x, y, subscripts, z, lower, upper, show.se, ...){
          if (grid) ticksGrid(x=1:length(levs), y=tickmarks$at)
          for (i in 1:length(zvals)){
            sub <- z[subscripts] == zvals[i]
            good <- !is.na(y[sub])
            os <- if(show.se)
              (i - (length(zvals) + 1)/2) * (2/(length(zvals)-1)) *
              .01 * (length(zvals) - 1) else 0
            effect.llines(x[sub][good]+os, y[sub][good], lwd=lwd, type='b', col=colors[.modc(i)],
                          pch=symbols[.mods(i)], lty=lines[.modl(i)], cex=cex, ...)
            if (show.se){
              larrows(x0=x[sub][good]+os, y0=lower[subscripts][sub][good],
                      x1=x[sub][good]+os, y1=upper[subscripts][sub][good],
                      angle=90, code=3, col=eval(colors[.modc(i)]),
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
        xlab=if (is.na(xlab)) predictors[x.var] else xlab,
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
      result <- update(plot, layout = if (is.null(layout))
        c(0, prod(dim(plot))) else layout)
      result$split <- split
      result$more <- more
      class(result) <- c("plot.eff", class(result))
    }
    ### multiline variate
    else{
      if (ci.style == "auto") ci.style <- "bands"
      effect.llines <- if(use.splines) spline.llines else effect.llines
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
      else range.adj(x[nm])
      tickmarks.x <- if ((nm %in% names(transform.x)) && !(is.null(transform.x))){
        trans <- transform.x[[nm]]$trans
        make.ticks(trans(xlm), link=transform.x[[nm]]$trans, inverse=transform.x[[nm]]$inverse, at=at, n=n)
      }
      else {
        trans <- I
        make.ticks(xlm, link=I, inverse=I, at=at, n=n)
      }
      key <- list(title=predictors[z.var], cex.title=1, border=TRUE,
                  text=list(as.character(zvals)),
                  lines=list(col=colors[.modc(1:length(zvals))], lty=lines[.modl(1:length(zvals))], lwd=lwd),
                  columns = if ("x" %in% names(key.args)) 1 else find.legend.columns(length(zvals)))
      for (k in names(key.args)) key[k] <- key.args[k]
      if (show.strip.values && n.predictors > 2){
        for (pred in predictors[-c(x.var, z.var)]){
          x[[pred]] <- as.factor(x[[pred]])
        }
      }
      plot <- xyplot(eval(parse(
        text=paste("fit ~trans(", predictors[x.var], ")",
                   if (n.predictors > 2) paste(" |",
                                               paste(predictors[-c(x.var, z.var)], collapse="*"))))),
        strip=function(...) strip.default(..., strip.names=c(factor.names, TRUE), sep=" = "),
        panel=function(x, y, subscripts, x.vals, rug, z, lower, upper, show.se, ...){
          if (grid) ticksGrid(x=tickmarks.x$at, y=tickmarks$at)
          if (rug && is.null(residuals)) lrug(trans(x.vals))
          axis.length <- diff(range(x))
          for (i in 1:length(zvals)){
            sub <- z[subscripts] == zvals[i]
            good <- !is.na(y[sub])
            effect.llines(x[sub][good], y[sub][good], lwd=lwd, type='l',
                          col=colors[.modc(i)], lty=lines[.modl(i)], cex=cex, ...)
            if(show.se){
              if(ci.style == "bars"){
                os <- (i - (length(zvals) + 1)/2) * (2/(length(zvals)-1)) *
                  .01 * axis.length
                larrows(x0=x[sub][good]+os, y0=lower[subscripts][sub][good],
                        x1=x[sub][good]+os, y1=upper[subscripts][sub][good],
                        angle=90, code=3, col=eval(colors[.modc(i)]),
                        length=.125*cex/1.5)
              }
              if(ci.style == "bands"){
                panel.bands(x[sub][good], y[sub][good],
                            upper[subscripts][sub][good], lower[subscripts][sub][good],
                            fill=eval(band.colors[.modb(i)]),
                            alpha=band.transparency, use.splines=use.splines)
              }
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
        xlab=if (is.na(xlab)) predictors[x.var] else xlab,
        x.vals=x.vals, rug=rug,
        z=x[,z.var],
        zvals=zvals,
        main=main,
        key=key,
        #
        lower=x$lower, upper=x$upper,
        show.se=has.se && ci.style %in% c("bars", "bands"),
        #
        data=x, scales=list(y=list(at=tickmarks$at, labels=tickmarks$labels),
                            rot=roty, x=list(at=tickmarks.x$at, labels=tickmarks.x$labels, rot=rotx),
                            alternating=alternating),  ...)
      result <- update(plot, layout = if (is.null(layout)) c(0, prod(dim(plot)))
                       else layout)
      result$split <- split
      result$more <- more
      class(result) <- c("plot.eff", class(result))
    }
    return(result)
  }
  # multiplot
  ci.style <- if(is.null(ci.style) || ci.style == "auto"){
    if(is.factor(x[, x.var])) "bars" else "bands"} else ci.style
  range <- if (has.se && ci.style !="none")
    range(c(x$lower, x$upper), na.rm=TRUE) else range(x$fit, na.rm=TRUE)
  # multiplot factor
  if (is.factor(x[,x.var])){
    
    ylim <- if (!any(is.na(ylim))) ylim else c(range[1] - .025*(range[2] - range[1]),
                                               range[2] + .025*(range[2] - range[1]))
    tickmarks <- if (type == "response" && rescale.axis) make.ticks(ylim, link=trans.link,
                                                                    inverse=trans.inverse, at=ticks$at, n=ticks$n)
    else make.ticks(ylim, link=I, inverse=I, at=ticks$at, n=ticks$n)
    
    levs <- levels(x[,x.var])
    if (show.strip.values){
      for (pred in predictors[-x.var]){
        x[[pred]] <- as.factor(x[[pred]])
      }
    }
    plot <- xyplot(eval(parse(
      text=paste("fit ~ as.numeric(", predictors[x.var], ") |",
                 paste(predictors[-x.var], collapse="*")))),
      strip=function(...) strip.default(..., strip.names=c(factor.names, TRUE), sep=" = "),
      panel=function(x, y, subscripts, lower, upper, has.se, ...){
        if (grid) ticksGrid(x=1:length(levs), y=tickmarks$at)
        good <- !is.na(y)
        if (has.se){
          if (ci.style == "bars"){
            larrows(x0=x[good], y0=lower[subscripts][good], x1=x[good], y1=upper[subscripts][good],
                    angle=90, code=3, col=colors[.modc(2)], length=0.125*cex/1.5)
          }
          else if(ci.style == "lines") {
            effect.llines(x[good], lower[subscripts][good], lty=2, col=colors[.modc(2)])
            effect.llines(x[good], upper[subscripts][good], lty=2, col=colors[.modc(2)])
          }
          else{ if(ci.style == "bands") {
            panel.bands(x[good], y[good], upper[subscripts][good], lower[subscripts][good],
                        fill=band.colors[1], alpha=band.transparency, use.splines=FALSE)
          }}
        }
        effect.llines(x[good], y[good], lwd=lwd, type='b', col=colors[1], pch=19, cex=cex, ...)
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
      xlab=if (is.na(xlab)) predictors[x.var] else xlab,
      scales=list(x=list(at=1:length(levs), labels=levs, rot=rotx),
                  y=list(at=tickmarks$at, labels=tickmarks$labels, rot=roty),
                  alternating=alternating),
      main=main,
      lower=x$lower, upper=x$upper, has.se=has.se, data=x, ...)
    result <- update(plot, layout = if (is.null(layout)) c(0, prod(dim(plot))) else layout)
    result$split <- split
    result$more <- more
    class(result) <- c("plot.eff", class(result))
  }
  ### multiplot variate   ***
  else{
    effect.llines <- if(use.splines) spline.llines else effect.llines
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
    else range.adj(x[nm])
    tickmarks.x <- if ((nm %in% names(transform.x)) && !(is.null(transform.x))){
      trans <- transform.x[[nm]]$trans
      make.ticks(trans(xlm), link=transform.x[[nm]]$trans, inverse=transform.x[[nm]]$inverse, at=at, n=n)
    }
    else {
      trans <- I
      make.ticks(xlm, link=I, inverse=I, at=at, n=n)
    }
    ylim <- if (!any(is.na(ylim))) ylim
    else if (is.null(residuals)) c(range[1] - .025*(range[2] - range[1]), range[2] + .025*(range[2] - range[1]))
    else if (rescale.axis) c(min(partial.residuals.range[1], range[1] - .025*(range[2] - range[1])),
                             max(partial.residuals.range[2], range[2] + .025*(range[2] - range[1])))
    else c(min(original.inverse(partial.residuals.range[1]), range[1] - .025*(range[2] - range[1])),
           max(original.inverse(partial.residuals.range[2]), range[2] + .025*(range[2] - range[1])))
    tickmarks <- if (type == "response" && rescale.axis) make.ticks(ylim, link=trans.link,
                                                                    inverse=trans.inverse, at=ticks$at, n=ticks$n)
    else make.ticks(ylim, link=I, inverse=I, at=ticks$at, n=ticks$n)
    x.fit <- x.data[, predictors[x.var]]
    use <- rep(TRUE, length(residuals))
    xx <- x[, predictors[-x.var], drop=FALSE]
    if (show.strip.values){
      for (pred in predictors[-x.var]){
        x[[pred]] <- as.factor(x[[pred]])
      }
    }
    plot <- xyplot(eval(parse(
      text=paste("fit ~ trans(", predictors[x.var], ") |",
                 paste(predictors[-x.var], collapse="*")))),
      strip=function(...) strip.default(..., strip.names=c(factor.names, TRUE), sep=" = "),
      panel=function(x, y, subscripts, x.vals, rug, lower, upper, has.se, ...){
        if (grid) ticksGrid(x=tickmarks.x$at, y=tickmarks$at)
        good <- !is.na(y)
        effect.llines(x[good], y[good], lwd=lwd, col=colors[1], ...)
        if (rug && is.null(residuals)) lrug(trans(x.vals))
        if (has.se){
          if (ci.style == "bars"){
            larrows(x0=x[good], y0=lower[subscripts][good],
                    x1=x[good], y1=upper[subscripts][good],
                    angle=90, code=3, col=eval(colors[.modc(2)]),
                    length=.125*cex/1.5)
          }
          else if(ci.style == "lines") {
            effect.llines(x[good], lower[subscripts][good], lty=2, col=colors[.modc(2)])
            effect.llines(x[good], upper[subscripts][good], lty=2, col=colors[.modc(2)])
          }
          else if(ci.style == "bands") {
            panel.bands(x[good], y[good], upper[subscripts][good], lower[subscripts][good],
                        fill=band.colors[1], alpha=band.transparency, use.splines=use.splines)
          }
        }
        if (!is.null(residuals)){
          predictors <- predictors[-x.var]
          factors <- sapply(xx, is.factor)
          for (predictor in predictors){
            use <- use & if(factors[predictor]) x.all[, predictor] == xx[subscripts[1], predictor]
            else x.all[, predictor] == xx[subscripts[1], predictor]
          }
          n.in.panel <- sum(use)
          if (n.in.panel > 0){
            fitted <- y[good][closest(trans(x.fit[use]), x[good])]
            partial.res <- if (!rescale.axis) original.inverse(original.link(fitted) + residuals[use])
            else fitted + residuals[use]
            lpoints(trans(x.fit[use]), partial.res, col=residuals.color, pch=residuals.pch, cex=residuals.cex)
            if (show.fitted) lpoints(trans(x.fit[use]), fitted, pch=16, col=residuals.color)  # REMOVE ME
            if (smooth.residuals && n.in.panel >= 10) {
              llines(loess.smooth(x.fit[use], partial.res, span=span, family=loess.family),
                     lwd=residuals.lwd, lty=residuals.lty, col=residuals.smooth.color)
            }
            if (id.n > 0){
              M <- cbind(trans(x.fit[use]), partial.res)
              md <- mahalanobis(M, colMeans(M), cov(M))
              biggest <- order(md, decreasing=TRUE)[1:id.n]
              pos <- ifelse(trans(x.fit[use][biggest]) > mean(current.panel.limits()$xlim), 2, 4)
              ltext(trans(x.fit[use][biggest]), partial.res[biggest], 
                    names(partial.res)[biggest], pos=pos, col=id.col, cex=id.cex)
            }
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
      xlab=if (is.na(xlab)) predictors[x.var] else xlab,
      x.vals=x.vals, rug=rug,
      main=main,
      lower=x$lower, upper=x$upper, has.se=has.se, data=x,
      scales=list(y=list(at=tickmarks$at, labels=tickmarks$labels, rot=roty),
                  x=list(at=tickmarks.x$at, labels=tickmarks.x$labels, rot=rotx),
                  alternating=alternating), ...)
    result <- update(plot, layout = if (is.null(layout)) c(0, prod(dim(plot))) else layout)
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

plot.efflist <- function(x, selection, rows, cols, ask=FALSE, graphics=TRUE, lattice, ...){
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
        print(plot(x[[(i-1)*cols + j]], lattice=lattice, ...))
      }
    }
  }
}
