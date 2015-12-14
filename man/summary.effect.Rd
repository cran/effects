\name{summary.eff}
\alias{print.eff}
\alias{print.effpoly}
\alias{print.efflatent}
\alias{print.efflist}
\alias{print.mlm.efflist}
\alias{print.summary.eff}
\alias{summary.eff}
\alias{summary.effpoly}
\alias{summary.efflatent}
\alias{summary.efflist}
\alias{summary.mlm.efflist}
\alias{plot.eff}
\alias{plot.effpoly}
\alias{plot.efflist}
\alias{plot.mlm.efflist}
\alias{setStrip}
\alias{restoreStrip}
\alias{[.efflist}

\title{Summarizing, Printing, and Plotting Effects}

\description{
\code{summary}, \code{print}, \code{plot}, and \code{[} methods for \code{eff}, \code{effpoly},
\code{efflist}, and \code{mlm.efflist} objects.
}

\usage{
\method{print}{eff}(x, type=c("response", "link"), ...)
\method{print}{effpoly}(x, type=c("probability", "logits"), ...)
\method{print}{efflatent}(x, ...)
\method{print}{efflist}(x, ...)
\method{print}{mlm.efflist}(x, ...)
\method{print}{summary.eff}(x, ...)
\method{summary}{eff}(object, type=c("response", "link"), ...)
\method{summary}{effpoly}(object, type=c("probability", "logits"), ...)
\method{summary}{efflatent}(object, ...)
\method{summary}{efflist}(object, ...)
\method{summary}{mlm.efflist}(object, ...)
\method{plot}{eff}(x, x.var,
    z.var=which.min(levels), multiline=is.null(x$se), rug=TRUE, 
    xlab, ylab, main=paste(effect, "effect plot"),
    colors=palette(), symbols=1:length(colors), lines=1:length(colors),
    cex=1.5, lwd=2, ylim, xlim=NULL,
    factor.names=TRUE, ci.style, band.transparency=0.15, band.colors=colors,
    type=c("rescale", "response", "link"), ticks=list(at=NULL, n=5),  
    alternating=TRUE, rotx=0, roty=0, grid=FALSE, layout, rescale.axis, 
    transform.x=NULL, ticks.x=NULL,
    key.args=NULL, 
    row=1, col=1, nrow=1, ncol=1, more=FALSE, 
    use.splines=TRUE, partial.residuals=c("adjusted", "raw"), show.fitted=FALSE,
    residuals.color="blue", residuals.pch=1, residuals.cex=1,
    smooth.residuals=TRUE, residuals.smooth.color=residuals.color, span=2/3, ...)
\method{plot}{effpoly}(x,
		type=c("probability", "logit"),
		x.var=which.max(levels),
		rug=TRUE,
		xlab,
		ylab=paste(x$response, " (", type, ")", sep=""),
		main=paste(effect, "effect plot"),
		colors, symbols, lines, cex=1.5, lwd=2,
		factor.names=TRUE, ci.style, band.colors, band.transparency=0.3,
		style=c("lines", "stacked"),
		confint=(style == "lines" && !is.null(x$confidence.level)),
		transform.x=NULL, ticks.x=NULL, xlim=NULL,
		ylim, rotx=0, alternating=TRUE, roty=0, grid=FALSE,
		layout, key.args=NULL,
		row=1, col=1, nrow=1, ncol=1, more=FALSE, use.splines=TRUE, ...)
\method{plot}{efflist}(x, selection, rows, cols, ask=FALSE, graphics=TRUE, ...)
\method{plot}{mlm.efflist}(x, ...)
setStrip(bg=3, fg="black", clip=c("off", "on"))
restoreStrip(saved)
\method{[}{efflist}(x, ...)
}

\arguments{
  \item{x}{an object of class \code{"eff"}, \code{"effpoly"}, \code{"efflist"}, \code{"mlm.efflist"},
  	or \code{"summary.eff"}, as appropriate.}
  \item{object}{an object of class \code{"eff"}, \code{"effpoly"}, \code{"efflist"}, or \code{"mlm.efflist"},
    as appropriate.}
  \item{type}{for printing or summarizing linear and generalized linear models, 
  	if \code{"response"} (the default), effects are printed on the scale of the response
    variable; if \code{"link"}, effects are printed  on the scale of the linear predictor.
    
    For plotting linear or genealized linearized models, \code{"rescale"} (the default) plots the vertical
    axis on the link scale (e.g., the logit scale for a logit model) but labels the axis on the response
    scale (e.g., the probability scale for a logit model); 
    \code{"response"} plots and labels the vertical axis on the scale of the response (e.g., the
    probability scale for a logit model); and
    \code{"link"} plots and labels the vertical axis
    on the scale of the link (e.g., the logit scale for a logit model).

    For polytomous logit models, this argument takes either \code{"probability"}
    or \code{"logit"}, with the former as the default.}
  \item{rescale.axis}{this argument is deprecated --- use the \code{type} argument instead. Setting 
    \code{rescale.axis=TRUE} is equivalent to \code{type="rescale"}; setting \code{rescale.axis=FALSE}
    is equivalent to \code{type="response"}. If specified, \code{rescale.axis} will override \code{type}.}
  \item{x.var}{the index (number) or quoted name of the covariate or factor to place on the 
    horizontal axis of each panel of the effect plot. The default is the 
    predictor with the largest number of levels or values.}
  \item{z.var}{for linear, generalized linear or mixed models,
  	the index (number) or quoted name of the covariate or factor for which
    individual lines are to be drawn in each panel of the effect plot. The default is the 
    predictor with the smallest number of levels or values. This argument is only
    used if \code{multiline = TRUE}.}
  \item{multiline}{for linear, generalized linear or mixed models,
  	if \code{TRUE}, each panel of the display represents combinations
    of values of two predictors, with one predictor (corresponding to \code{x.var})
    on the horzontal axis, and the other (corresponding to \code{z.var}) used to define
    lines in the graph; defaults to \code{TRUE} if there are no standard errors in
    the object being plotted, and \code{FALSE} otherwise.}
  \item{confint}{plot point-wise confidence bands around fitted effects (for
  	multinomial and proportional-odds logit models); defaults to \code{TRUE},
  	in which case separate panels are used for different response levels.}
  \item{rug}{if \code{TRUE}, the default, a rug plot is shown giving the \emph{marginal}
    distribution of the predictor on the horizontal axis, if this predictor is
    a covariate. The rug plot is suppressed if partial residuals are plotted.}
  \item{xlab}{the label for the horizontal axis of the effect plot; if missing, the
    function will use the name of the predictor on the horizontal axis.}
  \item{ylab}{the label for the vertical axis of the effect plot; the default is 
	constructed from the name of the
    response variable for the model from which the effect was computed.}
  \item{main}{the title for the plot, printed at the top; the default title is constructed from the
    name of the effect.}
  \item{colors}{\code{colors[1]} is used to plot effects, \code{colors[2]} to plot 
    confidence limits when \code{ci.style} is not equal to \code{"bands"}. In a mulitline plot, the successive \code{colors} correspond
    to the levels of the \code{z.var} covariate or factor. In a stacked plot or a plot
    without confidence bands for a
    multinomial or proportional-odds logit model, the successive \code{colors}
    correspond to the levels of the response factor. In all but stacked plots,
    \code{colors} defaults to \code{palette()}.   If \code{colors} is
    shorter than the number of levels, then it is recycled, so \code{colors="black"}
    will use black for all levels.

    For stacked multinomial-logit plots,
    \code{colors} defaults to \code{rainbow_hcl(levels)}, where \code{levels} is
    the number of levels of the response variable; for stacked proportional-odds model
    plots, \code{colors} defautls to \code{sequential_hcl(levels)}.   \code{colors}
    does not recycle for stacked plots.

    \bold{Warning:} This argument
    \emph{cannot} be abbreviated to \code{col}, which is used for a different
    purpose (see below).}
  \item{symbols, lines}{corresponding to the levels of the \code{z.var} covariate
    or factor on a multiline plot, or to the successive levels of the response factor
    in a line plot for a polytomous logit model. 
    These arguments are used only if \code{multiline = TRUE}
    or for polytomous logit models where the effects are plotted without confidence bands;
    in these cases a legend is drawn at the top of the display.  If these arguments are
    too short they are recycled.}
  \item{cex}{character expansion for plotted symbols; default is \code{1.5}.}
  \item{lwd}{line width for fitted lines.}
  \item{ylim}{2-element vector containing the lower and upper limits of the vertical axes;
    if \code{NULL}, the default, then the vertical axes are scaled from the data.}
  \item{xlim}{a named list of 2-element vectors, with the names corresponding to numeric
    predictors; if a numeric predictor is in the list, then when it appears on the horizontal
    axis, the axis limits will be taken from the corresponding vector; if a predictor is
    not in the list, or if the argument is \code{NULL} (the default), then the horizontal
    axis limits are computed from the data.}
  \item{factor.names}{a logical value, default \code{TRUE}, that controls the inclusion of
    factor names in conditioning-variable labels.}
  \item{ci.style}{confidence bounds can be indicated using error bars, using lines or confidence bands,
    depending on the plot type.
    For single line plots the default is \code{"bars"} for factors and \code{"bands"}
    for variates.  Style \code{"lines"} can also be used for either of these.  
    For multiline plots, the default is \code{"none"} for no 
    confidence bounds, but style \code{"bars"} or \code{"bands"} can also be used.  For a variate
    the option \code{"bars"} displays the error bars at
    each of the \code{levels} points for which the horizontal variable was evaluated.}
  \item{band.colors}{A vector of colors for the color of the confidence band with \code{ci.style="bands"}.
    The default is \code{band.colors=colors}.  For plots with one line, the choice, setting
    \code{band.colors="red"} produces a pleasing result, even if color provides no additional
    information in this case.  If this argument is too short it is recycled.}
  \item{band.transparency}{For \code{ci.style="bands"}, the alpha transparency of the fill color. Not all
    graphic devices support transparency.}
  \item{style}{(for multinomial or proportional-odds logit models) \code{"lines"} (the default
    for a line plot, or \code{"stacked"} for a stacked-bar or stacked-area plot. In the latter
    case only fitted probabilities may be plotted and confidence envelopes cannot be shown.}
  \item{ticks}{a two-item list controlling the placement of tick marks on the vertical axis,
    with elements \code{at} and \code{n}. If \code{at=NULL} (the default), the program
    attempts to find `nice' locations for the ticks, and the value of \code{n} (default, 
    \code{5}) gives the
    approximate number of tick marks desired; if \code{at} is non-\code{NULL}, then the
    value of \code{n} is ignored.}
  \item{ticks.x}{a named list of two-item lists controlling the placement of tick marks
    on the horizontal axis. Each list element is named for a numeric predictor in the model,
    and each sublist has elements \code{at} or \code{n} are for the \code{ticks} argument.
    If a predictor doesn't appear in the list, or if \code{ticks.x} is \code{NULL} (the default),
    then the tick marks are computed by the function.}
  \item{transform.x}{transformations to be applied to the horizontal axis, in the form of a named list,
    each of whose elements is itself a list of two functions, with sublist element names \code{trans} and \code{inverse}.
    The names of the list elements are numeric predictors in the model whose axes are to be transformed; the
    \code{trans} function is applied to the values of the predictor, and \code{inverse} is used for computing
    proper axis tick labels. If a numeric predictor is missing from \code{transform.x} then its axis is not
    transformed; if the argument is \code{NULL} (the default), then no predictor axes are transformed.}
  \item{alternating}{if \code{TRUE} (the default), the tick labels alternate by panels in
    multi-panel displays from left to right and top to bottom; if \code{FALSE}, tick labels
    appear at the bottom and on the left.}
  \item{rotx, roty}{rotation angles for the horizontal and vertical tick marks, 
    respectively.  Default is 0.}
  \item{grid}{if \code{TRUE}, add grid lines to the plot.  Default is FALSE.}
  \item{layout}{the \code{layout} argument to the \pkg{lattice} function \code{\link{xyplot}}
  	(or, in some cases \code{\link{densityplot}}), which
  	is used to draw the effect display; if not specified, the plot will be formatted so that
  	it appears on a single page.}
  \item{key.args}{additional arguments to be passed to the \code{key} trellis argument to 
    \code{\link{xyplot}} or \code{\link{densityplot}}, 
    e.g., to position the key (legend) in the plotting region; may also be used to modify
    the default behavior of the key.}
  \item{row, col, nrow, ncol, more}{These arguments are used to graph an effect as part of an
    array of plots; \code{row}, \code{col}, \code{nrow}, and \code{ncol} are used to compose
    the \code{split} argument and \code{more} the \code{more} argument to \code{\link{print.trellis}}.
    Normally these arguments are not set by the user, but by \code{plot.efflist}.
    \bold{Warning:} Note that \code{col} is \emph{not} used to specify colors;
    use \code{colors} instead (see above).}
  \item{selection}{the optional index (number) or quoted name of the effect in an effect
    list to be plotted; if not supplied, a menu of high-order terms is presented or all effects
    are plotted.}
  \item{rows, cols}{Number of rows and columns in the ``meta-array'' of plots produced for an \code{efflist} object;
  	if either argument is missing, then the meta-layout will be computed by the \code{plot} method.}
  \item{ask}{if \code{selection} is not supplied and \code{ask} is \code{TRUE},
    a menu of high-order terms is presented; if \code{ask} is \code{FALSE} (the default), effects for all
    high-order terms are plotted in an array.}
  \item{graphics}{if \code{TRUE} (the default), then the menu of terms to plot is presented
    in a dialog box rather than as a text menu.}
  \item{use.splines}{If \code{TRUE}, the default, then any lines drawn when the horizontal
    axis is not a factor use interpolating splines computed by the \code{\link{spline}} 
    function.  If FALSE, then linear interpoliation is used.  This argument is ignored if
    the horizontal axis is a factor.}
  \item{partial.residuals}{use \code{"adjusted"} partial residuals, computed at the panel-rounded
    values of the focal predictors not on the horizontal axis of each panel, or \code{"raw"}
    partial residuals, computed at the actual values of all focal predictors; the default is
    \code{"adjusted"}, and this argument is effective only when partial residuals are included
    in the effect object to be plotted --- see \code{\link{Effect}}.}
  \item{smooth.residuals}{whether to show a \code{\link{loess}} smooth of the partial residuals, if
    they are present; the default is \code{TRUE}. For a non-Gaussian \code{glm} model,
    a non-robust \code{loess} smooth is used; for a \code{lm} model or a Gaussian \code{glm} model,
    a robust smooth is employed.}
  \item{span}{of the \code{\link{xyplot}} smoother to be applied to partial residuals;
    the default is \code{2/3}.}
  \item{show.fitted}{if partial residuals are present in the effect object, also plot the partial
    fitted values (which will be shown as filled circles).}
  \item{residuals.color}{color for plotting partial residuals (default \code{"blue"}).}
  \item{residuals.smooth.color}{color for plotting the smooth of the partial residuals; the default is the \code{residuals.color}.}
  \item{residuals.pch}{plotting symbol for partial residuals (default \code{1}, open circles).}
  \item{residuals.cex}{character expansion (relative size) for symbols plotting partial residuals (default is \code{1}).}
  \item{bg}{if a single numeric value (the default is \code{3}), the color of the strips at the tops
    of lattice panels are set to gray scale, with the number of graditions, if there is more than one
    conditioning variable, corresponding to the value given (which will be rounded to a whole number).
    This argument can also be a vector of colors, specified in any manner recognized in R (e.g, by name).}
  \item{fg}{foreground color or colors for the strips at the top of lattice panels (the default is
    \code{"black"}); can be a single value or a vector of values of the same length as \code{bg}.}
  \item{clip}{\code{"off"} or \code{"on"}, determines whether or not conditioning values in the strips
    at the top of lattice panels are clipped on the left and right. The normal lattice default is
    \code{"on"}; the default in \code{setStrip} is \code{"off"}, allowing the lines representing numeric
    conditioning values to be displayed more clearly at the extreme left and right of strips.}
  \item{saved}{a set of lattice strip specifications returned by \code{setStrip}.}
  \item{...}{arguments to be passed down.}
}

\details{ 
  In a generalized linear model, by default, the \code{print} and \code{summary} methods for 
  \code{eff} objects print the computed effects on the scale of the
  response variable using the inverse of the 
  link function. In a logit model, for example, this means that the effects are expressed on the probability
  scale. 
  
  By default, effects in a GLM are plotted on the scale of the linear predictor, but the vertical
  axis is labelled on the response scale. This preserves the linear structure of the model while permitting
  interpretation on what is usually a more familiar scale.
  This approach may also be used with linear models, for example to display effects on the scale of the
  response even if the data are analyzed on a transformed scale, such as log or square-root.
  
  When a factor is on the x-axis, the \code{plot} method for \code{eff} objects
  connects the points representing the effect by line segments, creating a
  response ``profile.'' If you wish to suppress these lines, add the argument
  \code{lty=0} to the call to \code{plot} (see the examples).
  
  In a polytomous (multinomial or proportional-odds) logit model, by default effects are plotted on the
  probability scale; they may be alternatively plotted on the scale of the individual-level logits.
  
  The \code{setStrip} and \code{restoreStrip} functions modify the strips that appear in subsequent 
  lattice plots, including those produced by functions in the \pkg{effects} package. The default call
  \code{setStrip()} provides monochrome (rather than the lattice-default colored) strips with up to 3
  gray-scale values corresponding to 3 conditioning variables; clipping at the left and right of strips
  is also turned off by default by \code{setStrip}. \code{restoreStrip} may be used to reset
  lattice strips to previously saved parameters returned by \code{setStrip}.
}

\value{
  The \code{summary} method for \code{"eff"} objects returns a \code{"summary.eff"} object with the following components
  (those pertaining to confidence limits need not be present):
  \item{header}{a character string to label the effect.}
  \item{effect}{an array containing the estimated effect.}
  \item{lower.header}{a character string to label the lower confidence limits.}
  \item{lower}{an array containing the lower confidence limits.}
  \item{upper.header}{a character string to label the upper confidence limits.}
  \item{upper}{an array containing the upper confidence limits.}
  
  The \code{setStrip} function invisibly returns a list that can supply the argument of the 
\code{restoreStrip} function to restore the previous lattice strip specification.
  
  The \code{[} method for \code{"efflist"} objects is used to subset an \code{"efflist"} object and returns an object of the same class.
}
 
    
\author{John Fox \email{jfox@mcmaster.ca} and Jangman Hong.}


\seealso{\code{\link{effect}}, \code{\link{allEffects}}, \code{\link{xyplot}},
	\code{\link{densityplot}}, \code{\link{print.trellis}}, \code{\link{loess}},
	\code{\link[colorspace]{rainbow_hcl}}, \code{\link[colorspace:rainbow_hcl]{sequential_hcl}}}

\examples{
# also see examples in ?effect

mod.cowles <- glm(volunteer ~ sex + neuroticism*extraversion,
    data=Cowles, family=binomial)
eff.cowles <- allEffects(mod.cowles, xlevels=list(extraversion=seq(0, 24, 6)))
eff.cowles
as.data.frame(eff.cowles[[2]]) # neuroticism*extraversion interaction

plot(eff.cowles, 'sex', ylab="Prob(Volunteer)", grid=TRUE, rotx=90, lty=0)

.save.strip <- setStrip()

plot(eff.cowles, 'neuroticism:extraversion', ylab="Prob(Volunteer)",
    ticks=list(at=c(.1,.25,.5,.75,.9)))

\donttest{
# change color of the confidence bands to 'black' with .15 transparency
plot(eff.cowles, 'neuroticism:extraversion', ylab="Prob(Volunteer)",
    ticks=list(at=c(.1,.25,.5,.75,.9)), band.colors="red", band.transparency=.3)

plot(eff.cowles, 'neuroticism:extraversion', multiline=TRUE,
     ylab="Prob(Volunteer)", key.args = list(x = 0.75, y = 0.75, corner = c(0, 0)))

# use probability scale in place of logit scale, all lines are black.
plot(eff.cowles, 'neuroticism:extraversion', multiline=TRUE,
     ylab="Prob(Volunteer)", key.args = list(x = 0.75, y = 0.75, corner = c(0, 0)),
     colors="black", lines=1:8,
     ci.style="bands", type="response", band.colors=palette())

plot(effect('sex:neuroticism:extraversion', mod.cowles,
    xlevels=list(extraversion=seq(0, 24, 6))), multiline=TRUE)

plot(effect('sex:neuroticism:extraversion', mod.cowles,
    xlevels=list(extraversion=seq(0, 24, 6))), multiline=TRUE,
    type="response", ci.style="bands")
}

if (require(nnet)){
mod.beps <- multinom(vote ~ age + gender + economic.cond.national +
				economic.cond.household + Blair + Hague + Kennedy +
				Europe*political.knowledge, data=BEPS)
\donttest{
plot(effect("Europe*political.knowledge", mod.beps,
	xlevels=list(political.knowledge=0:3)))
}

plot(effect("Europe*political.knowledge", mod.beps,
		xlevels=list(political.knowledge=0:3),
		given.values=c(gendermale=0.5)),
	style="stacked", colors=c("blue", "red", "orange"), rug=FALSE)
}

if (require(MASS)){
mod.wvs <- polr(poverty ~ gender + religion + degree + country*poly(age,3),
	data=WVS)
plot(effect("country*poly(age, 3)", mod.wvs))

\donttest{
plot(effect("country*poly(age, 3)", mod.wvs), style="stacked",
	colors=c("gray75", "gray50", "gray25"))

plot(effect("country*poly(age, 3)", latent=TRUE, mod.wvs))
}
}

mod.pres <- lm(prestige ~ log(income, 10) + poly(education, 3) + poly(women, 2),
    data=Prestige)
eff.pres <- allEffects(mod.pres, default.levels=50)
\donttest{
plot(eff.pres)
plot(eff.pres[1:2])
}
plot(eff.pres[1],
    transform.x=list(income=list(trans=log10, inverse=function(x) 10^x)),
    ticks.x=list(income=list(at=c(1000, 2000, 5000, 10000, 20000))))
restoreStrip(.save.strip)
remove(.save.strip)
}

\keyword{hplot}
\keyword{models}
