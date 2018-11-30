# Plot method for effpoly objects

# modified by Michael Friendly: added ci.style="bands" & alpha.band= arg
# modified by Michael Friendly: added lwd= argument for llines (was lwd=2)
# 2013-11-06: fixed drop dimension when only one focal predictor. John
# 2014-10-10: namespace fixes. John
# 2014-12-05: made key.args more flexible. John
# 2014-03-22: use wide columns by default only when x for legend not set. J. Fox
# 2016-09-08: added show.strip.values argument to plot.effpoly(). J. Fox
# 2017-08-16: modified plot.effpoly() to consolidate arguments and use lattice theme. J. Fox
# 2017-08-20: reintroduce legacy arguments for plot.effpoly()
# 2017-08-20: introduced multiline argument under lines argument and as a "legacy" argument
# 2017-09-10: use replacement for grid.panel()
# 2017-11-22: added a check for non-estimable factor combinations with style="stacked"
# 2018-01-02, 2018-01-30: changed defaults for key.args, lines 140-141
# 2018-02-09: Use one-column key for stacked plot.
# 2018-02-28: Fix handling of rug arg (error reported by Dave Armstrong).
# 2018-07-08: add cex sub-args for x and y axes (suggestion of Charles Leger).
# 2018-07-08: add cex sub-arg for strips.
# 2018-10-05: modified plot.effpoly() so that multiline plots don't show confidence limits 
#             by default, and so that confidence bars for a factor are staggered.

plot.effpoly <- function(x, x.var=which.max(levels), 
                         main=paste(effect, "effect plot"),
                         symbols=TRUE, lines=TRUE, axes, confint, lattice, ...,
                         # legacy arguments:
                         type, multiline, rug, xlab, ylab, colors, cex, lty, lwd, 
                         factor.names, show.strip.values,
                         ci.style, band.colors, band.transparency, style, 
                         transform.x, ticks.x, xlim,
                         ticks, ylim, rotx, roty, alternating, grid, layout, 
                         key.args, use.splines){
  
  if (!is.logical(lines) && !is.list(lines)) lines <- list(lty=lines)
  lines <- applyDefaults(lines, 
              defaults=list(lty=trellis.par.get("superpose.line")$lty, 
                            lwd=trellis.par.get("superpose.line")$lwd[1], 
                            col=NULL, splines=TRUE, multiline=FALSE), 
              arg="lines")
  if (missing(multiline)) multiline <- lines$multiline
  if (missing(lwd)) lwd <- lines$lwd
  if (missing(use.splines)) use.splines <- lines$splines
  lines.col <- lines$col
  lines <- if (missing(lty)) lines$lty else lty
 
  if (!is.logical(symbols) && !is.list(symbols)) symbols <- list(pch=symbols)
  symbols <- applyDefaults(symbols,
                           defaults= list(
                             pch=trellis.par.get("superpose.symbol")$pch, 
                             cex=trellis.par.get("superpose.symbol")$cex[1]),
                           arg="symbols")
  cex <- symbols$cex
  symbols <- symbols$pch
  
  if (missing(axes)) axes <- NULL
  axes <- applyDefaults(axes, defaults=list(
    x=list(rotate=0, cex=1, rug=TRUE),
    y=list(lab=NULL, lim=c(NA, NA), ticks=list(at=NULL, n=5), 
           type="probability", rotate=0, cex=1),
    alternating=TRUE, grid=FALSE),
    arg="axes")
  
  x.args <- applyDefaults(axes$x, defaults=list(rotate=0, cex=1, rug=TRUE), 
                          arg="axes$x")
  
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
  cex.x <- x.args$cex
  x.args$rotate <- NULL
  x.args$rug <- NULL
  x.args$cex <- NULL
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
  if (length(xlab) == 0) xlab <- NULL
  if (length(xlim) == 0) xlim <- NULL
  if (length(ticks.x) == 0) ticks.x <- NULL
  if (length(transform.x) == 0) transform.x <- NULL
  
  y.args <- applyDefaults(axes$y, defaults=list(lab=NULL, lim=c(NA, NA), ticks=list(at=NULL, n=5), type="probability", style="lines", rotate=0, cex=1), arg="axes$y")
  if (missing(ylim)) ylim <- y.args$lim
  if (missing(ticks)) ticks <- y.args$ticks
  if (missing(type)) type <- y.args$type
  type <- match.arg(type, c("probability", "logit"))
  if (missing(ylab)) ylab <- y.args$lab
  if (is.null(ylab)) ylab <- paste0(x$response, " (", type, ")")
  if (missing(roty)) roty <- y.args$rotate
  cex.y <- y.args$cex
  if (missing(alternating)) alternating <- axes$alternating
  if (missing(grid)) grid <- axes$grid
  if (missing(style)) style <- match.arg(y.args$style, c("lines", "stacked"))
  
  if (missing(colors)) colors <- if (is.null(lines.col)){
    if (style == "lines" || x$model == "multinom")
      trellis.par.get("superpose.line")$col
    else sequential_hcl(length(x$y.levels))
  } else {
    lines.col
  }
  
  if (missing(confint)) confint <- NULL
  confint <- applyDefaults(confint,
                           defaults=list(style=if (style == "lines" && !multiline && !is.null(x$se.prob)) "auto" else "none", alpha=0.15, col=colors),
                           onFALSE=list(style="none", alpha=0, col="white"),
                           arg="confint")
  if (missing(ci.style)) ci.style <- confint$style
  if (missing(band.transparency)) band.transparency <- confint$alpha
  if (missing(band.colors)) band.colors <- confint$col
  if(!is.null(ci.style)) ci.style <- match.arg(ci.style, c("auto", "bars", "lines", "bands", "none")) 
  confint <- confint$style != "none"
  
  if (is.null(multiline)) multiline <- if (confint) FALSE else TRUE
  
  effect.llines <- llines
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
    ylim <- c(0, 1)
  }
  
  if (missing(lattice)) lattice <- NULL
  lattice <- applyDefaults(lattice, defaults=list(
    layout=NULL, #key.args=list(),  #New default added 1/2/2017 by sw
    strip=list(factor.names=TRUE, values=TRUE, cex=1),
    array=list(row=1, col=1, nrow=1, ncol=1, more=FALSE),
    arg="lattice"
  ))
  lattice$key.args <- applyDefaults(lattice$key.args, defaults=list(
    space="top", border=FALSE, fontfamily="sans", cex=.75, cex.title=1, 
    arg="key.args"
  ))
  if (missing(layout)) layout <- lattice$layout
  if (missing(key.args)) key.args <- lattice$key.args
  strip.args <- applyDefaults(lattice$strip, defaults=list(factor.names=TRUE, values=TRUE, cex=1), arg="lattice$strip")
  factor.names <- strip.args$factor.names
  if (missing(show.strip.values)) show.strip.values <- strip.args$values
  cex.strip <- strip.args$cex
  height.strip <- max(1, cex.strip)
  array.args <- applyDefaults(lattice$array, defaults=list(row=1, col=1, nrow=1, ncol=1, more=FALSE), arg="lattice$array")
  row <- array.args$row
  col <- array.args$col
  nrow <- array.args$nrow
  ncol <- array.args$ncol
  more <- array.args$more
  
  .mod <- function(a, b) ifelse( (d <- a %% b) == 0, b, d)
  .modc <- function(a) .mod(a, length(colors))
  .mods <- function(a) .mod(a, length(symbols))
  .modl <- function(a) .mod(a, length(lines))
  effect <- paste(sapply(x$variables, "[[", "name"), collapse="*")
  split <- c(col, row, ncol, nrow)
  n.predictors <- length(names(x$x))
  y.lev <- x$y.lev
  n.y.lev <- length(y.lev)
  ylevel.names <- make.names(paste("prob",y.lev))
  colnames(x$prob) <- colnames(x$logit) <- ylevel.names
  if (has.se){
    colnames(x$lower.logit) <- colnames(x$upper.logit) <- 
      colnames(x$lower.prob) <- colnames(x$upper.prob)<- ylevel.names
  }
  x.frame <-as.data.frame(x)
  predictors <- names(x.frame)[1:n.predictors]
  levels <- if (n.predictors==1) length (x.frame[,predictors])
  else sapply(apply(x.frame[, predictors, drop=FALSE], 2, unique), length)
  if (is.character(x.var)) {
    which.x <- which(x.var == predictors)
    if (length(which.x) == 0) stop(paste("x.var = '", x.var, "' is not in the effect.", sep=""))
    x.var <- which.x
  }
  x.vals <- x.frame[, names(x.frame)[x.var]]    
  response <- matrix(0, nrow=nrow(x.frame), ncol=n.y.lev)
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
  Data <- data.frame(prob, logit)
  if (has.se) Data <- cbind(Data, data.frame(lower.prob, upper.prob, lower.logit, upper.logit))
  Data[[x$response]] <- response
  for (i in 1:length(predictors)){
    Data <- cbind(Data, x.frame[predictors[i]])
  }
  levs <- levels(x$data[[predictors[x.var]]])
  n.predictor.cats <- sapply(Data[, predictors[-c(x.var)], drop=FALSE], 
                             function(x) length(unique(x)))
  if (length(n.predictor.cats) == 0) n.predictor.cats <- 1
  ci.style <- if(is.null(ci.style) || ci.style == "auto") {
    if(is.factor(x$data[[predictors[x.var]]])) "bars" else "bands"} else ci.style
  if( ci.style=="none" ) confint <- FALSE
  ### no confidence intervals if confint == FALSE or ci.style=="none"
  if (!confint){ # plot without confidence bands
    if (style == "lines"){ # line plot
      if (!multiline){
        layout <- if(is.null(layout)) c(prod(n.predictor.cats), length(levels(response)), 1) 
        else layout
        ### factor
        if (is.factor(x$data[[predictors[x.var]]])){ # x-variable a factor
          range <- if (type=="probability") range(prob, na.rm=TRUE) else range(logit, na.rm=TRUE)
          ylim <- if (!any(is.na(ylim))) ylim else c(range[1] - .025*(range[2] - range[1]),
                                                     range[2] + .025*(range[2] - range[1]))
          tickmarks <- make.ticks(ylim, link=I, inverse=I, at=ticks$at, n=ticks$n)
          levs <- levels(x$data[[predictors[x.var]]])
          if (show.strip.values){
            for (pred in predictors[-x.var]){
              Data[[pred]] <- as.factor(Data[[pred]])
            }
          }
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
            strip=strip.custom(strip.names=c(factor.names, TRUE), sep=" = ", par.strip.text=list(cex=cex.strip),
                               par.strip.text=list(cex=cex.strip)),
            panel=function(x, y, subscripts, x.vals, rug, ... ){
              if (grid) ticksGrid(x=1:length(levs), y=tickmarks$at)
              good <- !is.na(y)
              effect.llines(x[good], y[good], lwd=lwd, type="b", pch=19, col=colors[1], cex=cex, ...)
              subs <- subscripts+as.numeric(rownames(Data)[1])-1		
            },
            
            
            ylab=ylab,
            ylim=if (is.null(ylim))
              if (type == "probability") range(prob) else range(logit)
            else ylim,
            xlab=if (is.null(xlab)) predictors[x.var] else xlab[[x.var]],
            main=main,
            x.vals=x$data[[predictors[x.var]]],
            rug=rug,
            scales=list(x=list(at=1:length(levs), labels=levs, rot=rotx, cex=cex.x), 
                        y=list(at=tickmarks$at, labels=tickmarks$labels, rot=roty, cex=cex.y), alternating=alternating),
            layout=layout,
            data=Data, ...)
          result$split <- split
          result$more <- more
          class(result) <- c("plot.eff", class(result))
        }
        else { # x-variable numeric
          if(use.splines) effect.llines <- spline.llines # added 10/17/13
          range <- if (type=="probability") range(prob, na.rm=TRUE) else range(logit, na.rm=TRUE)
          ylim <- if (!any(is.na(ylim))) ylim else c(range[1] - .025*(range[2] - range[1]),
                                                     range[2] + .025*(range[2] - range[1]))
          tickmarks <- make.ticks(ylim, link=I, inverse=I, at=ticks$at, n=ticks$n)
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
          else range.adj(Data[nm]) # range(x.vals)
          tickmarks.x <- if ((nm %in% names(transform.x)) && !(is.null(transform.x))){
            trans <- transform.x[[nm]]$trans
            make.ticks(trans(xlm), link=transform.x[[nm]]$trans, inverse=transform.x[[nm]]$inverse, at=at, n=n)
          }
          else {
            trans <- I
            make.ticks(xlm, link=I, inverse=I, at=at, n=n)
          }
          if (show.strip.values){
            for (pred in predictors[-x.var]){
              Data[[pred]] <- as.factor(Data[[pred]])
            }
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
          strip=strip.custom(strip.names=c(factor.names, TRUE), sep=" = ", par.strip.text=list(cex=cex.strip),
                             par.strip.text=list(cex=cex.strip)),
          panel=function(x, y, subscripts, x.vals, rug, ... ){
            if (grid) ticksGrid(x=tickmarks.x$at, y=tickmarks$at)
            if (rug) lrug(trans(x.vals))
            good <- !is.na(y)
            effect.llines(x[good], y[good], lwd=lwd, col=colors[1], ...)
            subs <- subscripts+as.numeric(rownames(Data)[1])-1	
          },
          ylab=ylab,
          xlim=suppressWarnings(trans(xlm)),
          ylim= if (is.null(ylim))
            if (type == "probability") range(prob) else range(logit)
          else ylim,
          xlab=if (is.null(xlab)) predictors[x.var] else xlab[[x.var]],
          main=main,
          x.vals=x$data[[predictors[x.var]]],
          rug=rug,
          scales=list(y=list(at=tickmarks$at, labels=tickmarks$labels, rot=roty, cex=cex.y), 
                      x=list(at=tickmarks.x$at, labels=tickmarks.x$labels, rot=rotx, cex=cex.x),
                      alternating=alternating),
          layout=layout,
          data=Data, ...)
          result$split <- split
          result$more <- more
          class(result) <- c("plot.eff", class(result))
        }
        
      } else {
        layout <- if (is.null(layout)){
          lay <- c(prod(n.predictor.cats[-(n.predictors - 1)]), 
                   prod(n.predictor.cats[(n.predictors - 1)]), 1)
          if (lay[1] > 1) lay else lay[c(2, 1, 3)]
        }
        else layout
        if (n.y.lev > min(c(length(colors), length(lines), length(symbols))))
          warning('Colors, lines and symbols may have been recycled')
        range <- if (type=="probability") range(prob, na.rm=TRUE) else range(logit, na.rm=TRUE)
        ylim <- if (!any(is.na(ylim))) ylim else c(range[1] - .025*(range[2] - range[1]),
                                                   range[2] + .025*(range[2] - range[1]))
        tickmarks <- make.ticks(ylim, link=I, inverse=I, at=ticks$at, n=ticks$n)
        if (is.factor(x$data[[predictors[x.var]]])){ # x-variable a factor
          key <- list(title=x$response, cex.title=1, border=TRUE,
                      text=list(as.character(unique(response))),
                      lines=list(col=colors[.modc(1:n.y.lev)], lty=lines[.modl(1:n.y.lev)], lwd=lwd),
                      points=list(pch=symbols[.mods(1:n.y.lev)], col=colors[.modc(1:n.y.lev)]),
                      columns = if ("x" %in% names(key.args)) 1 else 
                        find.legend.columns(length(n.y.lev), 
                              space=if("x" %in% names(key.args)) "top" else key.args$space))
          for (k in names(key.args)) key[k] <- key.args[k]
          if (show.strip.values){
            for (pred in predictors[-x.var]){
              Data[[pred]] <- as.factor(Data[[pred]])
            }
          }
          result <- xyplot(eval(if (type=="probability") 
            parse(text=if (n.predictors==1) 
              paste("prob ~ as.numeric(", predictors[x.var], ")")
              else paste("prob ~ as.numeric(", predictors[x.var],") | ", 
                         paste(predictors[-x.var], collapse="*")))
            else parse(text=if (n.predictors==1) 
              paste("logit ~ as.numeric(", predictors[x.var], ")")
              else paste("logit ~ as.numeric(", predictors[x.var],") | ", 
                         paste(predictors[-x.var], collapse="*")))), 
            strip=strip.custom(strip.names=c(factor.names, TRUE), sep=" = ", par.strip.text=list(cex=cex.strip),
                               par.strip.text=list(cex=cex.strip)),
            panel=function(x, y, subscripts, rug, z, x.vals, ...){
              if (grid) ticksGrid(x=1:length(levs), y=tickmarks$at)
              for (i in 1:n.y.lev){
                sub <- z[subscripts] == y.lev[i]
                good <- !is.na(y[sub])
                effect.llines(x[sub][good], y[sub][good], lwd=lwd, type="b", col=colors[.modc(i)], lty=lines[.modl(i)],
                              pch=symbols[i], cex=cex, ...)
              }
            },
            ylab=ylab,
            ylim= if (is.null(ylim))
              if (type == "probability") range(prob) else range(logit)
            else ylim,
            xlab=if (is.null(xlab)) predictors[x.var] else xlab[[x.var]],
            x.vals=x$data[[predictors[x.var]]], 
            rug=rug,
            z=response,
            scales=list(x=list(at=1:length(levs), labels=levs, rot=rotx, cex=cex.x),
                        y=list(at=tickmarks$at, labels=tickmarks$labels, rot=roty, cex=cex.y),  
                        alternating=alternating),
            main=main,
            key=key,
            layout=layout,
            data=Data, ...)
          result$split <- split
          result$more <- more
          class(result) <- c("plot.eff", class(result))    		
        }
        else { # x-variable numeric
          if(use.splines) effect.llines <- spline.llines # added 10/17/13
          range <- if (type=="probability") range(prob, na.rm=TRUE) else range(logit, na.rm=TRUE)
          ylim <- if (!any(is.na(ylim))) ylim else c(range[1] - .025*(range[2] - range[1]),
                                                     range[2] + .025*(range[2] - range[1]))
          tickmarks <- make.ticks(ylim, link=I, inverse=I, at=ticks$at, n=ticks$n)
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
          else range.adj(Data[nm]) # range(x.vals)
          tickmarks.x <- if ((nm %in% names(transform.x)) && !(is.null(transform.x))){
            trans <- transform.x[[nm]]$trans
            make.ticks(trans(xlm), link=transform.x[[nm]]$trans, inverse=transform.x[[nm]]$inverse, at=at, n=n)
          }
          else {
            trans <- I
            make.ticks(xlm, link=I, inverse=I, at=at, n=n)
          }
          key <- list(title=x$response, cex.title=1, border=TRUE,
                      text=list(as.character(unique(response))), 
                      lines=list(col=colors[.modc(1:n.y.lev)], lty=lines[.modl(1:n.y.lev)], lwd=lwd),
                      columns = if ("x" %in% names(key.args)) 1 else 
                        find.legend.columns(length(n.y.lev), 
                                space=if("x" %in% names(key.args)) "top" else key.args$space))
          for (k in names(key.args)) key[k] <- key.args[k]
          if (show.strip.values){
            for (pred in predictors[-x.var]){
              Data[[pred]] <- as.factor(Data[[pred]])
            }
          }
          result <- xyplot(eval(if (type=="probability") 
            parse(text=if (n.predictors==1) paste("prob ~ trans(", predictors[x.var], ")")
                  else paste("prob ~ trans(", predictors[x.var],") |", 
                             paste(predictors[-x.var], collapse="*")))
            else parse(text=if (n.predictors==1) paste("logit ~ trans(", predictors[x.var], ")")
                       else paste("logit ~ trans(", predictors[x.var],") | ", 
                                  paste(predictors[-x.var], collapse="*")))), 
            strip=strip.custom(strip.names=c(factor.names, TRUE), sep=" = ", par.strip.text=list(cex=cex.strip),
                               par.strip.text=list(cex=cex.strip)),
            panel=function(x, y, subscripts, rug, z, x.vals, ...){
              if (grid) ticksGrid(x=tickmarks.x$at, y=tickmarks$at)
              if (rug) lrug(trans(x.vals))
              for (i in 1:n.y.lev){
                sub <- z[subscripts] == y.lev[i]
                good <- !is.na(y[sub])
                effect.llines(x[sub][good], y[sub][good], lwd=lwd, type="l", col=colors[.modc(i)], lty=lines[.modl(i)], ...)
              }
            },
            ylab=ylab,
            xlim=suppressWarnings(trans(xlm)),
            ylim= if (is.null(ylim))
              if (type == "probability") range(prob) else range(logit)
            else ylim,
            xlab=if (is.null(xlab)) predictors[x.var] else xlab[[x.var]],
            x.vals=x$data[[predictors[x.var]]], 
            rug=rug,
            z=response,
            scales=list(x=list(at=tickmarks.x$at, labels=tickmarks.x$labels, rot=rotx, cex=cex.x), 
                        y=list(at=tickmarks$at, labels=tickmarks$labels, rot=roty, cex=cex.y),
                        alternating=alternating),
            main=main,
            key=key,
            layout=layout,
            data=Data, ...)
          result$split <- split
          result$more <- more
          class(result) <- c("plot.eff", class(result))			
        }
      }
    }
    else { # stacked plot
      tickmarks <- make.ticks(c(0, 1), link=I, inverse=I, at=ticks$at, n=ticks$n)
      layout <- if (is.null(layout)){
        lay <- c(prod(n.predictor.cats[-(n.predictors - 1)]), 
                 prod(n.predictor.cats[(n.predictors - 1)]), 1)
        if (lay[1] > 1) lay else lay[c(2, 1, 3)]
      }
      else layout
      if (n.y.lev > length(colors))
        stop(paste('Not enough colors to plot', n.y.lev, 'regions'))
      key <- list(text=list(lab=rev(y.lev)), 
                  rectangle=list(col=rev(colors[1:n.y.lev])),
                  columns = 1)
                    # if ("x" %in% names(key.args)) 1 else 
                    # find.legend.columns(length(n.y.lev), 
                    #       space=if("x" %in% names(key.args)) "top" else key.args$space))
      for (k in names(key.args)) key[k] <- key.args[k]
      if (is.factor(x$data[[predictors[x.var]]])){ # x-variable a factor
# 11/22/17 check for rank deficient models and if found stop
        if(any(is.na(Data$prob))) stop("At least one combination of factor levels is not estimable.\n  Stacked plots are misleading, change to style='lines'")
        result <- barchart(eval(parse(text=if (n.predictors == 1) 
          paste("prob ~ ", predictors[x.var], sep="")
          else paste("prob ~ ", predictors[x.var]," | ", 
                     paste(predictors[-x.var], collapse="*")))), 
          strip=strip.custom(strip.names=c(factor.names, TRUE), sep=" = ", par.strip.text=list(cex=cex.strip),
                             par.strip.text=list(cex=cex.strip)),
          panel=function(x, y, ...){
            panel.barchart(x, y, ...)
            if (grid) ticksGrid(x=NA, y=tickmarks$at, col="white")
          },
          groups = response,
          col=colors,
          horizontal=FALSE, 
          stack=TRUE, 
          data=Data, 
          ylim=ylim, # if (is.null(ylim)) 0:1 else ylim,
          ylab=ylab, 
          xlab=if (is.null(xlab)) predictors[x.var] else xlab[[x.var]],
          scales=list(x=list(rot=rotx, at=1:length(levs), labels=levs, cex=cex.x), 
                      y=list(rot=roty, at=tickmarks$at, labels=tickmarks$labels, cex=cex.y), 
                      alternating=alternating),
          main=main,
          key=key,
          layout=layout)
        result$split <- split
        result$more <- more
        class(result) <- c("plot.eff", class(result))			
      }
      else { # x-variable numeric
        if(use.splines) effect.llines <- spline.llines # added 10/17/13
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
        else range.adj(Data[nm]) # range(x.vals)
        tickmarks.x <- if ((nm %in% names(transform.x)) && !(is.null(transform.x))){
          trans <- transform.x[[nm]]$trans
          make.ticks(trans(xlm), link=transform.x[[nm]]$trans, inverse=transform.x[[nm]]$inverse, at=at, n=n)
        }
        else {
          trans <- I
          make.ticks(xlm, link=I, inverse=I, at=at, n=n)
        }
        if (show.strip.values){
          for (pred in predictors[-x.var]){
            x$x[[pred]] <- as.factor(x$x[[pred]])
          }
        }
        result <- densityplot(eval(parse(text=if (n.predictors == 1)
          paste("~ trans(", predictors[x.var], ")", sep="")
          else paste("~ trans(", predictors[x.var], ") | ",
                     paste(predictors[-x.var], collapse="*")))),
          probs=x$prob,
          strip=strip.custom(strip.names=c(factor.names, TRUE), sep=" = ", par.strip.text=list(cex=cex.strip),
                             par.strip.text=list(cex=cex.strip)),
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
            if (rug) lrug(trans(x.vals))
            if (grid) ticksGrid(x=tickmarks.x$at, y=tickmarks$at, col="white")
          },
          rug=rug,
          x.vals=x$data[[predictors[x.var]]],
          data=x$x,
          xlim=suppressWarnings(trans(xlm)),
          ylim= c(0, 1), # if (is.null(ylim)) 0:1 else ylim,
          ylab=ylab,
          xlab=if (is.null(xlab)) predictors[x.var] else xlab[[x.var]],
          scales=list(x=list(at=tickmarks.x$at, labels=tickmarks.x$labels, rot=rotx, cex=cex.x), 
                      y=list(rot=roty, at=tickmarks$at, labels=tickmarks$labels, cex=cex.y),
                      alternating=alternating),
          main=main,
          key=key,
          layout=layout, ...)
        result$split <- split
        result$more <- more
        class(result) <- c("plot.eff", class(result))
      }
    }
  }
  ### with confidence bands
  else{ # plot with confidence bands
    if (type == "probability"){
      lower <- lower.prob
      upper <- upper.prob
    }
    else {
      lower <- lower.logit
      upper <- upper.logit
    }
    if (!multiline){
      layout <- if(is.null(layout)) c(prod(n.predictor.cats), length(levels(response)), 1) 
      else layout
      ### factor
      if (is.factor(x$data[[predictors[x.var]]])){ # x-variable a factor
        range <- range(c(lower, upper), na.rm=TRUE)
        ylim <- if (!any(is.na(ylim))) ylim else c(range[1] - .025*(range[2] - range[1]),
                                                   range[2] + .025*(range[2] - range[1]))
        tickmarks <- make.ticks(ylim, link=I, inverse=I, at=ticks$at, n=ticks$n)
        levs <- levels(x$data[[predictors[x.var]]])
        if (show.strip.values){
          for (pred in predictors[-x.var]){
            Data[[pred]] <- as.factor(Data[[pred]])
          }
        }
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
          strip=strip.custom(..., strip.names=c(factor.names, TRUE), sep=" = ", par.strip.text=list(cex=cex.strip),
                             par.strip.text=list(cex=cex.strip)),
          panel=function(x, y, subscripts, x.vals, rug, lower, upper, ... ){
            if (grid) ticksGrid(x=1:length(levs), y=tickmarks$at)
            good <- !is.na(y)
            effect.llines(x[good], y[good], lwd=lwd, type="b", pch=19, col=colors[1], cex=cex, ...)
            subs <- subscripts+as.numeric(rownames(Data)[1])-1		
            if (ci.style == "bars"){
              larrows(x0=x[good], y0=lower[subs][good], 
                      x1=x[good], y1=upper[subs][good], 
                      angle=90, code=3, col=colors[.modc(2)], length=0.125*cex/1.5)
            }
            else  if(ci.style == "lines"){
              effect.llines(x[good], lower[subs][good], lty=2, col=colors[.modc(2)])
              effect.llines(x[good], upper[subs][good], lty=2, col=colors[.modc(2)])
            }
            else { if(ci.style == "bands") {		
              panel.bands(x[good], y[good],
                          lower[subs][good], upper[subs][good],
                          fill=band.colors[1], alpha=band.transparency)
            }}
          },
          
          
          ylab=ylab,
          ylim= if (is.null(ylim)) c(min(lower), max(upper)) else ylim,
          xlab=if (is.null(xlab)) predictors[x.var] else xlab[[x.var]],
          main=main,
          x.vals=x$data[[predictors[x.var]]],
          rug=rug,
          lower=lower,
          upper=upper, 
          scales=list(x=list(at=1:length(levs), labels=levs, rot=rotx, cex=cex.x), 
                      y=list(at=tickmarks$at, labels=tickmarks$labels, rot=roty, cex=cex.y), 
                      alternating=alternating),
          layout=layout,
          data=Data, ...)
        result$split <- split
        result$more <- more
        class(result) <- c("plot.eff", class(result))
      }
      else { # x-variable numeric
        if(use.splines) effect.llines <- spline.llines # added 10/17/13
        range <- range(c(lower, upper), na.rm=TRUE)
        ylim <- if (!any(is.na(ylim))) ylim else c(range[1] - .025*(range[2] - range[1]),
                                                   range[2] + .025*(range[2] - range[1]))
        tickmarks <- make.ticks(ylim, link=I, inverse=I, at=ticks$at, n=ticks$n)
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
        else range.adj(Data[nm]) # range(x.vals)
        tickmarks.x <- if ((nm %in% names(transform.x)) && !(is.null(transform.x))){
          trans <- transform.x[[nm]]$trans
          make.ticks(trans(xlm), link=transform.x[[nm]]$trans, inverse=transform.x[[nm]]$inverse, at=at, n=n)
        }
        else {
          trans <- I
          make.ticks(xlm, link=I, inverse=I, at=at, n=n)
        }
        if (show.strip.values){
          for (pred in predictors[-x.var]){
            Data[[pred]] <- as.factor(Data[[pred]])
          }
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
        strip=strip.custom(strip.names=c(factor.names, TRUE), sep=" = ", par.strip.text=list(cex=cex.strip),
                           par.strip.text=list(cex=cex.strip)),
        panel=function(x, y, subscripts, x.vals, rug, lower, upper, ... ){
          if (grid) ticksGrid(x=tickmarks.x$at, y=tickmarks$at)
          if (rug) lrug(trans(x.vals))
          good <- !is.na(y)
          effect.llines(x[good], y[good], lwd=lwd, col=colors[1], ...)
          subs <- subscripts+as.numeric(rownames(Data)[1])-1	
          if (ci.style == "bars"){
            larrows(x0=x[good], y0=lower[subs][good], 
                    x1=x[good], y1=upper[subs][good], 
                    angle=90, code=3, col=colors[.modc(2)], length=0.125*cex/1.5)
          }
          else  if(ci.style == "lines"){
            effect.llines(x[good], lower[subs][good], lty=2, col=colors[.modc(2)])
            effect.llines(x[good], upper[subs][good], lty=2, col=colors[.modc(2)])
          } 
          else { if(ci.style == "bands") {	
            panel.bands(x[good], y[good],
                        lower[subs][good], upper[subs][good],
                        fill=band.colors[1], alpha=band.transparency)
          }}
        },
        ylab=ylab,
        xlim=suppressWarnings(trans(xlm)),
        ylim= if (is.null(ylim)) c(min(lower), max(upper)) else ylim,
        xlab=if (is.null(xlab)) predictors[x.var] else xlab[[x.var]],
        main=main,
        x.vals=x$data[[predictors[x.var]]],
        rug=rug,
        lower=lower,
        upper=upper, 
        scales=list(y=list(at=tickmarks$at, labels=tickmarks$labels, rot=roty, cex=cex.y), 
                    x=list(at=tickmarks.x$at, labels=tickmarks.x$labels, rot=rotx, cex=cex.x),
                    alternating=alternating),
        layout=layout,
        data=Data, ...)
        result$split <- split
        result$more <- more
        class(result) <- c("plot.eff", class(result))
      }
    } else {
      
      layout <- if (is.null(layout)){
        lay <- c(prod(n.predictor.cats[-(n.predictors - 1)]), 
                 prod(n.predictor.cats[(n.predictors - 1)]), 1)
        if (lay[1] > 1) lay else lay[c(2, 1, 3)]
      }
      else layout
      if (n.y.lev > min(c(length(colors), length(lines), length(symbols))))
        warning('Colors, lines and symbols may have been recycled')
      if (is.factor(x$data[[predictors[x.var]]])){ # x-variable a factor
        range <- range(c(lower, upper), na.rm=TRUE)
        ylim <- if (!any(is.na(ylim))) ylim else c(range[1] - .025*(range[2] - range[1]),
                                                   range[2] + .025*(range[2] - range[1]))
        tickmarks <- make.ticks(ylim, link=I, inverse=I, at=ticks$at, n=ticks$n)
        key <- list(title=x$response, cex.title=1, border=TRUE,
                    text=list(as.character(unique(response))),
                    lines=list(col=colors[.modc(1:n.y.lev)], lty=lines[.modl(1:n.y.lev)], lwd=lwd),
                    points=list(pch=symbols[.mods(1:n.y.lev)], col=colors[.modc(1:n.y.lev)]),
                    columns = if ("x" %in% names(key.args)) 1 else 
                      find.legend.columns(length(n.y.lev), 
                            space=if("x" %in% names(key.args)) "top" else key.args$space))

        for (k in names(key.args)) key[k] <- key.args[k]
        if (show.strip.values){
          for (pred in predictors[-x.var]){
            Data[[pred]] <- as.factor(Data[[pred]])
          }
        }
        result <- xyplot(eval(if (type=="probability") 
          parse(text=if (n.predictors==1) 
            paste("prob ~ as.numeric(", predictors[x.var], ")")
            else paste("prob ~ as.numeric(", predictors[x.var],") | ", 
                       paste(predictors[-x.var], collapse="*")))
          else parse(text=if (n.predictors==1) 
            paste("logit ~ as.numeric(", predictors[x.var], ")")
            else paste("logit ~ as.numeric(", predictors[x.var],") | ", 
                       paste(predictors[-x.var], collapse="*")))), 
          strip=strip.custom(strip.names=c(factor.names, TRUE), sep=" = ", par.strip.text=list(cex=cex.strip),
                             par.strip.text=list(cex=cex.strip)),
          panel=function(x, y, subscripts, rug, z, x.vals, lower, upper, ...){
            if (grid) ticksGrid(x=1:length(levs), y=tickmarks$at)
            for (i in 1:n.y.lev){
              os <- if (ci.style == "bars"){
                (i - (n.y.lev + 1)/2) * (2/(n.y.lev-1)) * .01 * (n.y.lev - 1)
              } else {
                0
              }
              sub <- z[subscripts] == y.lev[i]
              good <- !is.na(y[sub])
              effect.llines(x[sub][good] + os, y[sub][good], lwd=lwd, type="b", col=colors[.modc(i)], lty=lines[.modl(i)],
                            pch=symbols[i], cex=cex, ...)
              if (ci.style == "bars"){
                larrows(x0=x[sub][good] + os, y0=lower[ ][sub][good], 
                        x1=x[sub][good] + os, y1=upper[subscripts][sub][good], 
                        angle=90, code=3, col=colors[.modc(i)], length=0.125*cex/1.5)
              }
              else  if(ci.style == "lines"){
                effect.llines(x[sub][good], lower[subscripts][sub][good], lty=lines[.modl(i)], col=colors[.modc(i)])
                effect.llines(x[sub][good], upper[subscripts][sub][good], lty=lines[.modl(i)], col=colors[.modc(i)])
              }
              else { if(ci.style == "bands") {		
                panel.bands(x[sub][good], y[sub][good],
                            lower[subscripts][sub][good], upper[subscripts][sub][good],
                            fill=colors[.modc(i)], alpha=band.transparency)
              }}
            }
          },
          ylab=ylab,
          ylim= if (is.null(ylim)) c(min(lower), max(upper)) else ylim,
          xlab=if (is.null(xlab)) predictors[x.var] else xlab[[x.var]],
          x.vals=x$data[[predictors[x.var]]], 
          rug=rug,
          z=response,
          lower=lower,
          upper=upper,
          scales=list(x=list(at=1:length(levs), labels=levs, rot=rotx, cex=cex.x),
                      y=list(at=tickmarks$at, labels=tickmarks$labels, rot=roty, cex=cex.y),  
                      alternating=alternating),
          main=main,
          key=key,
          layout=layout,
          data=Data, ...)
        result$split <- split
        result$more <- more
        class(result) <- c("plot.eff", class(result))    		
      }
      else { # x-variable numeric
        if(use.splines) effect.llines <- spline.llines # added 10/17/13
        range <- range(c(lower, upper), na.rm=TRUE)
        ylim <- if (!any(is.na(ylim))) ylim else c(range[1] - .025*(range[2] - range[1]),
                                                   range[2] + .025*(range[2] - range[1]))
        tickmarks <- make.ticks(ylim, link=I, inverse=I, at=ticks$at, n=ticks$n)
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
        else range.adj(Data[nm]) # range(x.vals)
        tickmarks.x <- if ((nm %in% names(transform.x)) && !(is.null(transform.x))){
          trans <- transform.x[[nm]]$trans
          make.ticks(trans(xlm), link=transform.x[[nm]]$trans, inverse=transform.x[[nm]]$inverse, at=at, n=n)
        }
        else {
          trans <- I
          make.ticks(xlm, link=I, inverse=I, at=at, n=n)
        }
        key <- list(title=x$response, cex.title=1, border=TRUE,
                    text=list(as.character(unique(response))), 
                    lines=list(col=colors[.modc(1:n.y.lev)], lty=lines[.modl(1:n.y.lev)], lwd=lwd),
                    columns = if ("x" %in% names(key.args)) 1 else 
                      find.legend.columns(length(n.y.lev), 
                            space=if("x" %in% names(key.args)) "top" else key.args$space))
        for (k in names(key.args)) key[k] <- key.args[k]
        if (show.strip.values){
          for (pred in predictors[-x.var]){
            Data[[pred]] <- as.factor(Data[[pred]])
          }
        }
        result <- xyplot(eval(if (type=="probability") 
          parse(text=if (n.predictors==1) paste("prob ~ trans(", predictors[x.var], ")")
                else paste("prob ~ trans(", predictors[x.var],") |", 
                           paste(predictors[-x.var], collapse="*")))
          else parse(text=if (n.predictors==1) paste("logit ~ trans(", predictors[x.var], ")")
                     else paste("logit ~ trans(", predictors[x.var],") | ", 
                                paste(predictors[-x.var], collapse="*")))), 
          strip=strip.custom(strip.names=c(factor.names, TRUE), sep=" = ", par.strip.text=list(cex=cex.strip),
                             par.strip.text=list(cex=cex.strip)),
          panel=function(x, y, subscripts, rug, z, x.vals, lower, upper, ...){
            if (grid) ticksGrid(x=tickmarks.x$at, y=tickmarks$at)
            if (rug) lrug(trans(x.vals))
            for (i in 1:n.y.lev){
              sub <- z[subscripts] == y.lev[i]
              good <- !is.na(y[sub])
              effect.llines(x[sub][good], y[sub][good], lwd=lwd, type="l", col=colors[.modc(i)], lty=lines[.modl(i)], ...)
              if (ci.style == "bars"){
                larrows(x0=x[sub][good], y0=lower[subscripts][sub][good], 
                        x1=x[sub][good], y1=upper[subscripts][sub][good], 
                        angle=90, code=3, col=colors[.modc(i)], length=0.125*cex/1.5)
              }
              else  if(ci.style == "lines"){
                effect.llines(x[sub][good], lower[subscripts][sub][good], lty=lines[.modl(i)], col=colors[.modc(i)])
                effect.llines(x[sub][good], upper[subscripts][sub][good], lty=lines[.modl(i)], col=colors[.modc(i)])
              }
              else { if(ci.style == "bands") {		
                panel.bands(x[sub][good], y[sub][good],
                            lower[subscripts][sub][good], upper[subscripts][sub][good],
                            fill=colors[.modc(i)], alpha=band.transparency)
              }}
            }
          },
          ylab=ylab,
          xlim=suppressWarnings(trans(xlm)),
          ylim= if (is.null(ylim)) c(min(lower), max(upper)) else ylim,
          xlab=if (is.null(xlab)) predictors[x.var] else xlab[[x.var]],
          x.vals=x$data[[predictors[x.var]]], 
          rug=rug,
          z=response,
          lower=lower,
          upper=upper,
          scales=list(x=list(at=tickmarks.x$at, labels=tickmarks.x$labels, rot=rotx, cex=cex.x), 
                      y=list(at=tickmarks$at, labels=tickmarks$labels, rot=roty, cex=cex.y),
                      alternating=alternating),
          main=main,
          key=key,
          layout=layout,
          data=Data, ...)
        result$split <- split
        result$more <- more
        class(result) <- c("plot.eff", class(result))			
      }
      
    }
  }
  result
}
