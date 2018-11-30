\name{predictorEffects}
\alias{predictorEffect}
\alias{predictorEffect.poLCA}
\alias{predictorEffect.svyglm}
\alias{predictorEffect.default}
\alias{predictorEffects}
\alias{predictorEffects.poLCA}
\alias{predictorEffects.default}
\title{
Functions For Computing Predictor Effects
}
\description{
Alternatives to the \code{Effect} and \code{allEffects} functions that use a different paradigm for conditioning in an effect display.  The user specifies one predictor, either continuous or a factor, for the horizontal axis of a plot, and the function determines the appropriate plot to display (which is drawn by \code{plot}). See the vignette \href{../doc/predictor-effects-gallery.pdf}{Predictor Effects Graphics Gallery} for details and examples.
}
\usage{
predictorEffect(predictor, mod, focal.levels=50, xlevels=5, ...)

\method{predictorEffect}{poLCA}(predictor, mod, focal.levels=50, 
    xlevels=5, ...)

\method{predictorEffect}{svyglm}(predictor, mod, focal.levels=50, 
    xlevels=5, ...)

\method{predictorEffect}{default}(predictor, mod, focal.levels=50, 
    xlevels=5, ...)

predictorEffects(mod, predictors, focal.levels=50, xlevels=5, ...)

\method{predictorEffects}{poLCA}(mod, predictors = ~ ., 
    focal.levels=50, xlevels=5, ...)

\method{predictorEffects}{default}(mod, predictors = ~ .,
    focal.levels=50, xlevels=5,  ...)

}
\arguments{
  \item{mod}{
A model object.  Supported models include all those described on the help page for \code{\link{Effect}}.
}
  \item{predictor}{quoted name of the focal predictor.}
  \item{predictors}{
If the default \code{~ .}, a predictor effect plot is drawn for each predictor (not regressor) in a model.  Otherwise, this should be a one-sided formula listing the first-order predictors for which predictor effect plots should be drawn.
}
  \item{focal.levels}{for \code{predictorEffect}, the number of evenly-spaced values (the default is 50) for the focal predictor or a vector of values for the
  focal predictor. 
  
  For \code{predictorEffects}, the number of evenly-spaced values (default 50) to use for each focal predictor in turn or a named list,
  similar to \code{xlevels}, giving the number of values or the values themselves for each predictor individually, to be used when that predictor is the
  focal predictor; if a focal predictor doesn't appear in the list, the default of 50 values is used.}
  
  \item{xlevels}{this argument is used to set the levels of conditioning predictors; it may either
  be a single number specifying the number of evenly-spaced values (the default is 5) to which each conditioning predictor is to be set, 
  or it may be a list with elements named for the predictors giving the number of values or a vector of values to which each
  conditioning predictor is to be set, as explained in the help for \code{\link{Effect}}. If the
  focal predictor is included in the \code{xlevels} list, it is disregarded; if any conditioning predictor
  is omitted from the list, its number of values is set to 5. 
  
  The \code{xlevels} argument works similarly for \code{predictorEffect} and \code{predictorEffects}.}
  \item{\dots}{
Additional arguments passed to \code{\link{Effect}}.
}
}
\details{
Effect plots view a fitted regression function E(Y|X) in (sequences of) two-dimensional plots using conditioning and slicing.  The functions described here use a different method of determining the conditioning and slicing than \code{allEffects} uses.  The predictor effect of a focal predictor, say \code{x1},is the  usual effect for the generalized interaction of \code{x1} with all the other predictors in a model.  When a predictor effect object is plotted, the focal predictor is by default plotted on the horizontal axis.

For example, in the model \code{mod} with formula \code{y ~ x1 + x2 + x3}, the predictor effect \code{p1 <- predictorEffects(mod, ~ x1)} is essentially equilavent to \code{p2 <- Effect("x1", mod)}.  When plotted, these objects may produce different graphs because \code{plot(p1)} will always put \code{x1} on the horizontal axis while \code{plot(p2)} uses a rule to determine the horizontal axis based on the characteristics of all the predictors, e.g., preferring continuous predictors over factors.

If \code{mod} has the formula \code{y ~ x1 + x2 + x3 + x1:x2}, then \code{p1 <- predictorEffects(mod, ~ x1)} is essentially equivalent to \code{p2 <- Effect(c("x1", "x2"), mod)}.  As in the last example, the plotted versions of these objects may differ because of different rules used to determine the predictor on the horizontal axis.

If \code{mod} has the formula \code{y ~ x1 + x2 + x3 + x1:x2 + x1:x3}, then \code{p1 <- predictorEffects(mod, ~ x1)} is essentially equilavent to \code{p2 <- Effect(c("x1", "x2", "x3"), mod)}.  Again, the plotted versions of these objects may differ because of the rules used to determine the horizontal axis.
}
\value{
\code{predictorEffect} returns an object of class \code{c("predictoreff", "eff")}.  The components of the object are described in the help for \code{\link{Effect}}; \code{predictorEffects} returns an object of class \code{"predictorefflist"}, which is a list whose elements are of class \code{c("predictoreff", "eff")}.
}
\references{
See \code{\link{Effect}}.
}
\author{
S. Weisberg \email{sandy@umn.edu} and J. Fox
}



\seealso{
\code{\link{Effect}}, \code{\link{plot.predictoreff}}, the \href{../doc/predictor-effects-gallery.pdf}{Predictor Effects Graphics Gallery} vignette, and the \href{../doc/partial-residuals.pdf}{Effect Displays with Partial Residuals} vignette.
}
\examples{
mod <- lm(prestige ~ type*(education + income) + women, Prestige)
plot(predictorEffect("income", mod))
plot(predictorEffects(mod, ~ education + income + women))

mod.cowles <- glm(volunteer ~ sex + neuroticism*extraversion, data=Cowles, family=binomial)
plot(predictorEffects(mod.cowles, xlevels=4))
plot(predictorEffect("neuroticism", mod.cowles, xlevels=list(extraversion=seq(5, 20, by=5))),
     axes=list(grid=TRUE,
               x=list(rug=FALSE),
               y=list(lab="Probability of Vounteering")),
     lines=list(multiline=TRUE), 
     type="response")
predictorEffects(mod.cowles, focal.levels=4, xlevels=4)

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
