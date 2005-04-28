\name{summary.effect}
\alias{print.effect}
\alias{print.effect.list}
\alias{print.summary.effect}
\alias{summary.effect}
\alias{summary.effect.list}
\alias{plot.effect}
\alias{plot.effect.list}

\title{Summarizing, Printing, and Plotting Effects}

\description{
\code{summary}, \code{print}, and \code{plot} methods for \code{effect}
and \code{effect.list} objects.
}

\usage{
\method{print}{effect}(x, type=c("response", "link"), ...)
\method{print}{effect.list}(x, ...)
\method{print}{summary.effect}(x, ...)
\method{summary}{effect}(object, type=c("response", "link"), ...)
\method{summary}{effect.list}(object, ...)
\method{plot}{effect}(x, x.var=which.max(levels),
    z.var=which.min(levels), multiline=is.null(x$se), rug=TRUE, xlab,
    ylab=x$response, main=paste(effect, "effect plot"),
    colors=palette(), symbols=1:10, lines=1:10, cex=1.5, ylim,
    factor.names=TRUE, type=c("response", "link"), ticks=list(at=NULL, n=5), 
    alternating=TRUE, rescale.axis=TRUE, 
    row=1, col=1, nrow=1, ncol=1, more=FALSE, ...)
\method{plot}{effect.list}(x, selection, ask=TRUE, ...)
}

\arguments{
  \item{x}{an object of type \code{"effect"}, \code{"summary.effect"}, or \code{"effect.list"},
    as appropriate.}
  \item{object}{an object of type \code{"effect"} or \code{"effect.list"},
    as appropriate.}
  \item{type}{if \code{"response"} (the default), effects are printed 
    or the vertical axis is labelled on the scale of the response
    variable; if \code{"link"}, effects are printed 
    or the vertical axis labelled on the scale of the linear predictor.}
  \item{x.var}{the index (number) or quoted name of the covariate or factor to place on the 
    horizontal axis of each panel of the effect plot. The default is the 
    predictor with the largest number of levels or values.}
  \item{z.var}{the index (number) or quoted name of the covariate or factor for which
    individual lines are to be drawn in each panel of the effect plot. The default is the 
    predictor with the smallest number of levels or values. This argument is only
    used if \code{multiline = TRUE}.}
  \item{multiline}{if \code{TRUE}, each panel of the display represents combinations
    of values of two predictors, with one predictor (corresponding to \code{x.var})
    on the horzontal axis, and the other (corresponding to \code{z.var}) used to define
    lines in the graph; defaults to \code{TRUE} if there are no standard errors in
    the object being plotted, and \code{FALSE} otherwise.}
  \item{rug}{if \code{TRUE}, the default, a rug plot is shown giving the marginal 
    distribution of the predictor on the horizontal axis, if this predictor is
    a covariate.}
  \item{xlab}{the label for the horizontal axis of the effect plot; if missing, the
    function will use the name of the predictor on the horizontal axis.}
  \item{ylab}{the label for the vertical axis of the effect plot; the default is the
    response variable for the model from which the effect was computed.}
  \item{main}{the title for the plot, printed at the top; the default title is constructed from the
    name of the effect.}
  \item{colors}{\code{colors[1]} is used to plot effects, \code{colors[2]} to plot 
    confidence bands. In a mulitline plot, the successive \code{colors} correspond
    to the levels of the \code{z.var} covariate or factor.}
  \item{symbols, lines}{corresponding to the levels of the \code{z.var} covariate
    or factor on a multiline plot. These arguments are used only if \code{multiline = TRUE};
    in this case a legend is drawn at the top of the display.}
  \item{cex}{character expansion for plotted symbols; default is \code{1}.}
  \item{ylim}{2-element vector containing the lower and upper limits of the vertical axes;
    if \code{NULL}, the default, then the vertical axes are scaled from the data.}
  \item{factor.names}{a logical value, default \code{TRUE}, that controls the inclusion of
    factor names in conditioning-variable labels.}
  \item{ticks}{a two-item list controlling the placement of tick marks on the vertical axis,
    with elements \code{at} and \code{n}. If \code{at=NULL} (the default), the program
    attempts to find `nice' locations for the ticks, and the value of \code{n} (default, 
    \code{5}) gives the
    approximate number of tick marks desired; if \code{at} is non-\code{NULL}, then the
    value of \code{n} is ignored.}
  \item{alternating}{if \code{TRUE} (the default), the tick labels alternate by panels in
    multi-panel displays from left to right and top to bottom; if \code{FALSE}, tick labels
    appear at the bottom and on the left.}
  \item{rescale.axis}{if \code{TRUE} (the default), the tick marks on the vertical axis
    are labelled on the response scale (e.g., the probability scale for effects computed on
    the logit scale for a binomial GLM).}
  \item{row, col, nrow, ncol, more}{These arguments are used to graph an effect as part of an
    array of plots; \code{row}, \code{col}, \code{nrow}, and \code{ncol} are used to compose
    the \code{split} argument and \code{more} the \code{more} argument to \code{\link[lattice]{print.trellis}}.
    Normally these arguments are not set by the user, but by \code{print.effect.list}.}
  \item{selection}{the optional index (number) or quoted name of the effect in an effect
    list to be plotted; if not supplied, a menu of high-order terms is presented or all effects
    are plotted.}
  \item{ask}{if \code{selection} is not supplied and \code{ask} is \code{TRUE} (the default),
    a menu of high-order terms is presented; if \code{ask} is \code{FALSE}, effects for all
    high-order terms are plotted in an array.}
  \item{...}{arguments to be passed down.}
}

\details{ 
  In a generalized linear model, by default, the \code{print} and \code{summary} methods for 
  \code{effect} objects print the computed effects on the scale of the
  response variable using the inverse of the 
  link function. In a logit model, for example, this means that the effects are expressed on the probability
  scale. 
  
  By default, effects in a GLM are plotted on the scale of the linear predictor, but the vertical
  axis is labelled on the response scale. This preserves the linear structure of the model while permitting
  interpretation on what is usually a more familiar scale.
  This approach may also be used with linear models, for example to display effects on the scale of the
  response even if the data are analyzed on a transformed scale, such as log or square-root.
}

\value{
  The \code{summary} method for \code{"effect"} objects returns a \code{"summary.effect"} object with the following components
  (those pertaining to confidence limits need not be present):
  \item{header}{a character string to label the effect.}
  \item{effect}{an array containing the estimated effect.}
  \item{lower.header}{a character string to label the lower confidence limits.}
  \item{lower}{an array containing the lower confidence limits.}
  \item{upper.header}{a character string to label the upper confidence limits.}
  \item{upper}{an array containing the upper confidence limits.}
}
 
    
\author{John Fox \email{jfox@mcmaster.ca}.}


\seealso{\code{\link{effect}}, \code{\link{all.effects}}, \code{\link[lattice]{xyplot}}}

\examples{
data(Cowles)
mod.cowles <- glm(volunteer ~ sex + neuroticism*extraversion, 
    data=Cowles, family=binomial)
eff.cowles <- all.effects(mod.cowles, xlevels=list(neuroticism=0:24, 
    extraversion=seq(0, 24, 6)))
eff.cowles
    \dontrun{
    model: volunteer ~ sex + neuroticism * extraversion
    
    sex effect
    sex
    female      male 
    0.4409441 0.3811941 
    
    neuroticism*extraversion effect
            extraversion
    neuroticism          0         6        12        18        24
            0  0.07801066 0.1871263 0.3851143 0.6301824 0.8225756
            1  0.08636001 0.1963396 0.3870453 0.6200668 0.8083638
            2  0.09551039 0.2058918 0.3889798 0.6098458 0.7932997
            3  0.10551835 0.2157839 0.3909179 0.5995275 0.7773775
           . . .
            23 0.51953129 0.4747277 0.4303273 0.3870199 0.3454282
            24 0.54709527 0.4895731 0.4323256 0.3768303 0.3243880
    }
plot(eff.cowles, 'sex', ylab="Prob(Volunteer)")
    \dontrun{
    Loading required package: lattice 
    }

plot(eff.cowles, 'neuroticism:extraversion', ylab="Prob(Volunteer)",
    ticks=list(at=c(.1,.25,.5,.75,.9)))

plot(eff.cowles, 'neuroticism:extraversion', multiline=TRUE, 
    ylab="Prob(Volunteer)")
    
plot(effect('sex:neuroticism:extraversion', mod.cowles,
    xlevels=list(neuroticism=0:24, extraversion=seq(0, 24, 6))), multiline=TRUE)
\dontrun{
    Warning message: 
    sex:neuroticism:extraversion does not appear in the model in: 
      effect("sex:neuroticism:extraversion", mod.cowles, 
      xlevels = list(neuroticism = 0:24,
}

data(Prestige)
mod.pres <- lm(prestige ~ log(income, 10) + poly(education, 3) + poly(women, 2), 
    data=Prestige)
eff.pres <- all.effects(mod.pres, default.levels=50)
plot(eff.pres, ask=FALSE)

}

\keyword{hplot}
\keyword{models}
