\name{summary.eff}
\alias{print.eff}
\alias{print.effpoly}
\alias{print.efflatent}
\alias{print.efflist}
\alias{print.summary.eff}
\alias{summary.eff}
\alias{summary.effpoly}
\alias{summary.efflatent}
\alias{summary.efflist}
\alias{plot.eff}
\alias{plot.effpoly}
\alias{plot.efflist}
\alias{[.efflist}

\title{Summarizing, Printing, and Plotting Effects}

\description{
\code{summary}, \code{print}, \code{plot}, and \code{[} methods for \code{eff}, \code{effpoly},
and \code{efflist} objects.
}

\usage{
\method{print}{eff}(x, type=c("response", "link"), ...)
\method{print}{effpoly}(x, type=c("probability", "logits"), ...)
\method{print}{efflatent}(x, ...)
\method{print}{efflist}(x, ...)
\method{print}{summary.eff}(x, ...)
\method{summary}{eff}(object, type=c("response", "link"), ...)
\method{summary}{effpoly}(object, type=c("probability", "logits"), ...)
\method{summary}{efflatent}(object, ...)
\method{summary}{efflist}(object, ...)
\method{plot}{eff}(x, x.var = which.max(levels), z.var = which.min(levels), 
    multiline = is.null(x$se), rug = TRUE, xlab, ylab, main = paste(effect, 
        "effect plot"), colors = palette(), symbols = 1:10, lines = 1:10, 
    cex = 1.5, ylim, xlim = NULL, factor.names = TRUE, 
    ci.style, 
    type = c("response", "link"), ticks = list(at = NULL, n = 5), alternating = TRUE, 
    rotx = 0, roty = 0, grid = FALSE, layout, rescale.axis = TRUE, 
    transform.x = NULL, ticks.x = NULL, key.args = NULL, row = 1, 
    col = 1, nrow = 1, ncol = 1, more = FALSE, ...) 
\method{plot}{effpoly}(x, type = c("probability", "logit"), x.var = which.max(levels), 
    rug = TRUE, xlab, ylab = paste(x$response, " (", type, ")", 
        sep = ""), main = paste(effect, "effect plot"), colors, 
    symbols = 1:10, lines = 1:10, cex = 1.5, factor.names = TRUE,
    ci.style, 
    style = c("lines", "stacked"), confint = (style == "lines" && 
        !is.null(x$confidence.level)), transform.x = NULL, ticks.x = NULL, 
    xlim = NULL, ylim, rotx = 0, alternating = TRUE, roty = 0, 
    grid = FALSE, layout, key.args = NULL, row = 1, col = 1, 
    nrow = 1, ncol = 1, more = FALSE, ...)
\method{plot}{efflist}(x, selection, rows, cols, ask=FALSE, graphics=TRUE, ...)
\method{[}{efflist}(x, ...)
}

\arguments{
  \item{x}{an object of class \code{"eff"}, \code{"effpoly"}, \code{"efflist"}, 
  	or \code{"summary.eff"}, as appropriate.}
  \item{object}{an object of class \code{"eff"}, \code{"effpoly"}, or \code{"efflist"},
    as appropriate.}
  \item{type}{for linear and generalized linear models, 
  	if \code{"response"} (the default), effects are printed 
    or the vertical axis is labelled on the scale of the response
    variable; if \code{"link"}, effects are printed 
    or the vertical axis labelled on the scale of the linear predictor.
    For polytomous logit models, this argument takes either \code{"probability"}
    or \code{"logit"}, with the former as the default.}
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
  \item{rug}{if \code{TRUE}, the default, a rug plot is shown giving the marginal 
    distribution of the predictor on the horizontal axis, if this predictor is
    a covariate.}
  \item{xlab}{the label for the horizontal axis of the effect plot; if missing, the
    function will use the name of the predictor on the horizontal axis.}
  \item{ylab}{the label for the vertical axis of the effect plot; the default is 
	constructed from the name of the
    response variable for the model from which the effect was computed.}
  \item{main}{the title for the plot, printed at the top; the default title is constructed from the
    name of the effect.}
  \item{colors}{\code{colors[1]} is used to plot effects, \code{colors[2]} to plot 
    confidence bands. In a mulitline plot, the successive \code{colors} correspond
    to the levels of the \code{z.var} covariate or factor. In a stacked plot or a plot
    without confidence bands for a
    multinomial or proportional-odds logit model, the successive \code{colors}
    correspond to the levels of the response factor. In all but stacked plots,
    \code{colors} defaults to \code{palette()}; for stacked multinomial-logit plots,
    \code{colors} defaults to \code{rainbow_hcl(levels)}, where \code{levels} is
    the number of levels of the response variable; for stacked proportional-odds model
    plots, \code{colors} defautls to \code{sequential_hcl(levels)}.
    \bold{Warning:} This argument
    \emph{cannot} be abbreviated to \code{col}, which is used for a different
    purpose (see below).}
  \item{symbols, lines}{corresponding to the levels of the \code{z.var} covariate
    or factor on a multiline plot, or to the successive levels of the response factor
    in a line plot for a polytomous logit model. 
    These arguments are used only if \code{multiline = TRUE}
    or for polytomous logit models where the effects are plotted without confidence bands;
    in these cases a legend is drawn at the top of the display.}
  \item{cex}{character expansion for plotted symbols; default is \code{1.5}.}
  \item{ylim}{2-element vector containing the lower and upper limits of the vertical axes;
    if \code{NULL}, the default, then the vertical axes are scaled from the data.}
  \item{xlim}{a named list of 2-element vectors, with the names corresponding to numeric
    predictors; if a numeric predictor is in the list, then when it appears on the horizontal
    axis, the axis limits will be taken from the corresponding vector; if a predictor is
    not in the list, or if the argument is \code{NULL} (the default), then the horizontal
    axis limits are computed from the data.}
  \item{factor.names}{a logical value, default \code{TRUE}, that controls the inclusion of
    factor names in conditioning-variable labels.}
  \item{ci.style}{confidence bounds can be indicated using error bars or using lines.
    For single line plots the default is \code{"bars"} for factors and \code{"lines"}
    for variates.  For multiline plots, the default is \code{"none"} for no 
    confidence bounds, but style \code{"bars"} can also be used.  For a variate
    the option \code{"bars"} displays the error bars at
    each of the \code{levels} points for which the horizontal variable was evaluated.}
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
  \item{rescale.axis}{if \code{TRUE} (the default), the tick marks on the vertical axis
    are labelled on the response scale (e.g., the probability scale for effects computed on
    the logit scale for a binomial GLM).}
  \item{key.args}{additional arguments to be passed to the \code{key} trellis argument to 
    \code{\link{xyplot}} or \code{\link{densityplot}}, 
    e.g., to position the key (legend) in the plotting region.}
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
  
  In a polytomous (multinomial or proportional-odds) logit model, by default effects are plotted on the
  probability scale; they may be alternatively plotted on the scale of the individual-level logits.
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
  
  The \code{[} method for \code{"efflist"} objects is used to subset an \code{"efflist"} object and returns an object of the same class.
}
 
    
\author{John Fox \email{jfox@mcmaster.ca} and Jangman Hong.}


\seealso{\code{\link{effect}}, \code{\link{allEffects}}, \code{\link{xyplot}},
	\code{\link{densityplot}}, \code{\link{print.trellis}}
	\code{\link[colorspace]{rainbow_hcl}}, \code{\link[colorspace:rainbow_hcl]{sequential_hcl}}}

\examples{
mod.cowles <- glm(volunteer ~ sex + neuroticism*extraversion, 
    data=Cowles, family=binomial)
eff.cowles <- allEffects(mod.cowles, xlevels=list(neuroticism=0:24, 
    extraversion=seq(0, 24, 6)))
eff.cowles

plot(eff.cowles, 'sex', ylab="Prob(Volunteer)", grid=TRUE, rotx=90)

plot(eff.cowles, 'neuroticism:extraversion', ylab="Prob(Volunteer)",
    ticks=list(at=c(.1,.25,.5,.75,.9)))

plot(eff.cowles, 'neuroticism:extraversion', multiline=TRUE, 
    ylab="Prob(Volunteer)", key.args = list(x = 0.75, y = 0.75, corner = c(0, 0)))
    
plot(effect('sex:neuroticism:extraversion', mod.cowles,
    xlevels=list(neuroticism=0:24, extraversion=seq(0, 24, 6))), multiline=TRUE)

if (require(nnet)){
mod.beps <- multinom(vote ~ age + gender + economic.cond.national +
				economic.cond.household + Blair + Hague + Kennedy +
				Europe*political.knowledge, data=BEPS)
plot(effect("Europe*political.knowledge", mod.beps, 
	xlevels=list(Europe=1:11, political.knowledge=0:3)))
	
plot(effect("Europe*political.knowledge", mod.beps, 
		xlevels=list(Europe=1:11, political.knowledge=0:3),
		given.values=c(gendermale=0.5)),
	style="stacked", colors=c("blue", "red", "orange"), rug=FALSE)
}
	
if (require(MASS)){
mod.wvs <- polr(poverty ~ gender + religion + degree + country*poly(age,3),
	data=WVS)
plot(effect("country*poly(age, 3)", mod.wvs))

plot(effect("country*poly(age, 3)", mod.wvs), style="stacked", 
	colors=c("gray75", "gray50", "gray25"))
	
plot(effect("country*poly(age, 3)", latent=TRUE, mod.wvs))
}

mod.pres <- lm(prestige ~ log(income, 10) + poly(education, 3) + poly(women, 2), 
    data=Prestige)
eff.pres <- allEffects(mod.pres, default.levels=50)
plot(eff.pres)
plot(eff.pres[1:2])
plot(eff.pres[1], 
    transform.x=list(income=list(trans=log10, inverse=function(x) 10^x)),
    ticks.x=list(income=list(at=c(1000, 2000, 5000, 10000, 20000))))

}

\keyword{hplot}
\keyword{models}
