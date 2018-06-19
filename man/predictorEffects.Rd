\name{predictorEffects}
\alias{predictorEffect}
\alias{predictorEffect.svyglm}
\alias{predictorEffect.default}
\alias{predictorEffects}
\alias{predictorEffects.default}
\title{
Functions For Computing Predictor Effects
}
\description{
Alternatives to the \code{Effect} and \code{allEffects} functions that use a different paradigm for conditioning in an effects display.  The user specifies one predictor, either continuous or a factor, for the horizontal axis of a plot, and the function determines the appropriate plot to display (which is drawn by \code{plot}).
}
\usage{
predictorEffect(predictor, mod, xlevels, ...)

\method{predictorEffect}{svyglm}(predictor, mod, xlevels, ...)

\method{predictorEffect}{default}(predictor, mod, xlevels, ...)

predictorEffects(mod, predictors, ...)

\method{predictorEffects}{default}(mod, predictors = ~ ., ...)

}
\arguments{
  \item{mod}{
A model object.  Supported models include all those described on the help page for \code{\link{Effect}}.
}
  \item{predictor}{quoted name of the focal predictor.}
  \item{xlevels}{this argument is used to set the values for any predictor in the
  effect that is not a factor.  For a predictor
  effect, the default is to use 50 quantiles of the focal predictor on the x-axis
  between the 0.01 and 0.98 quantiles. See \code{\link{Effect}} for details and for
  how other defaults are set.}
  \item{predictors}{
If the default \code{~ .}, a predictor effecrts plot is drawn for each predictor (not regressor) in a model.  Otherwise, this should be a one-sided formula listing the first-order predictors for which predictor effects plots should be drawn.
}
  \item{\dots}{
Additional arguments passed to \code{\link{Effect}}.
}
}
\details{
Effects plots view a fitted regression function E(Y|X) in (sequences of) two-dimensional plots using conditioning and slicing.  The functions describe here use a different method of determining the conditioning and slicing than \code{Effects} uses.  The predictor effects a focal predictor say \code{x1} will be the the usual effect for the generalized interaction of \code{x1} with all the other predictors in a model.  When a predictor effects object is plotted, the focal predictor is by default plotted on the horizontal axis.

For example, in the model \code{mod} with formula \code{y ~ x1 + x2 + x3}, then \code{p1 <- predictorEffects(mod, ~ x1)} is essentially equilavent to \code{p2 <- Effect("x1", mod)}.  When plotted, these objects may be different because \code{plot(p1)} will always put \code{x1} on the horizontal axis while \code{plot(p2)} uses a rule to determine the horizontal axis based on the characteristics of all the predictors, preferring continuous predictors over factors.

If \code{mod} has the formula \code{y ~ x1 + x2 + x3 + x1:x2}, then \code{p1 <- predictorEffects(mod, ~ x1)} is essentially equilavent to \code{p2 <- Effect(c("x1", "x2"), mod)}.  As in the last example, the plotted versions of these objects may differ because of rules used to determine the horizontal axis.

If \code{mod} has the formula \code{y ~ x1 + x2 + x3 + x1:x2 + x1:x3}, then \code{p1 <- predictorEffects(mod, ~ x1)} is essentially equilavent to \code{p2 <- Effect(c("x1", "x2", "x3"), mod)}.  The plotted versions of these objects may differ because of rules used to determine the horizontal axis.
}
\value{
\code{predictorEffect} returns an object of class \code{c(predictoreff, eff)}.  The components of the object are described under the detalis at \code{\link{Effect}}.  \code{predictorEffects} returns an object of class \code{predictorefflist}, which is a list whose elements are of class \code{c(predictoreff, eff)}
}
\references{
See \code{\link{Effect}}.
}
\author{
S. Weisberg, \email{sandy@umn.edu}
}



\seealso{
\code{\link{Effect}}, \code{\link{plot.predictoreff}}
}
\examples{
mod <- lm(prestige ~ type*(education + income) + women, Prestige)
plot(predictorEffect("income", mod))
plot(predictorEffects(mod, ~ education + income + women))

# svyglm() example (adapting an example from the survey package)
\donttest{
if (require(survey)){
  data(api)
  dstrat<-svydesign(id=~1, strata=~stype, weights=~pw,
    data=apistrat, fpc=~fpc)
  mod <- svyglm(sch.wide ~ ell + meals + mobility, design=dstrat,
    family=quasibinomial())
  plot(predictorEffects(mod),
    axes=list(y=list(lim=log(c(0.4, 0.99)/c(0.6, 0.01)),
      ticks=list(at=c(0.4, 0.75, 0.9, 0.95, 0.99)))))
}
}
}

\keyword{hplot}
\keyword{models}
