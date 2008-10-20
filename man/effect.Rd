\name{effect}
\alias{effect}
\alias{effect.lm}
\alias{effect.multinom}
\alias{effect.polr}
\alias{allEffects}
\alias{as.data.frame.eff}
\alias{as.data.frame.effpoly}
\alias{as.data.frame.efflatent}

\title{Functions For Constructing Effect Plots}
\description{
  \code{effect} constructs an \code{"eff"} object for a term (usually a high-order term) 
  in a linear or generalized linear model, or an \code{"effpoly"} object for a term in a
   multinomial or proportional-odds logit model, 
  absorbing the lower-order terms marginal
  to the term in question, and averaging over other terms in the model. 
  
  \code{allEffects} identifies all of the high-order terms in a model and returns
  a list of \code{"eff"} or \code{"effpoly"} objects (i.e., an object of type \code{"efflist"}).
}

\usage{
effect(term, mod, ...)

\method{effect}{lm}(term, mod, xlevels=list(), default.levels=10, given.values, 
	se=TRUE, confidence.level=.95, 
    transformation=list(link=family(mod)$linkfun, inverse=family(mod)$linkinv), 
    typical=mean, ...)
    
\method{effect}{multinom}(term, mod, confidence.level=.95, xlevels=list(), default.levels=10, 
	given.values, se=TRUE, typical=mean, ...)
		
\method{effect}{polr}(term, mod, confidence.level=.95, xlevels=list(), default.levels=10, 
	given.values, se=TRUE, typical=mean, latent=FALSE, ...)
    
allEffects(mod, ...)

\method{as.data.frame}{eff}(x, row.names=NULL, optional=TRUE, ...)

\method{as.data.frame}{effpoly}(x, row.names=NULL, optional=TRUE, ...)

\method{as.data.frame}{efflatent}(x, row.names=NULL, optional=TRUE, ...)
}

\arguments{
  \item{term}{the quoted name of a term, usually, but not necessarily, a high-order 
    term in the model. The term must be given exactly as it appears in the printed
    model, although either colons (\code{:}) or asterisks (\code{*}) may be used
    for interactions.}
  \item{mod}{an object of class \code{"lm"}, \code{"glm"}, \code{"multinom"}, or \code{"polr"}.}
  \item{xlevels}{an optional list of values at which to set covariates,
    with components of the form \code{covariate.name = vector.of.values}.}
  \item{default.levels}{number of values for covariates that are not
    specified explicitly via \code{xlevels}; covariate values set by
    default are evenly spaced between the minimum and maximum values in
    the data.}
  \item{given.values}{a numeric vector of named elements, setting particular
  	columns of the model matrix to specific values for terms \emph{not}
  	appearing in an effect; if specified, takes precedence over the
  	application of the function given in the \code{typical} argument
  	(below). Care must be taken in specifying these values --- e.g.,
  	for a factor, the values of all contrasts should be given and these
  	should be consistent with each other.}
  \item{se}{if \code{TRUE}, the default, calculate standard errors and
    confidence limits for the effects.}
  \item{confidence.level}{level at which to compute confidence limits
    based on the standard-normal distribution; the default is \code{0.95}.}
  \item{transformation}{a two-element list with elements \code{link} and
    \code{inverse}. For a generalized linear model, these are by default
    the link function and inverse-link (mean) function. For a linear model,
    these default to \code{NULL}. If \code{NULL}, the identify function,
    \code{I}, is used; this effect can also be achieved by setting the
    argument to \code{NULL}. The inverse-link may be used to transform effects
    when they are printed or plotted; the link may be used in positioning
    axis labels (see below). If the link is not given, an attempt will be
    made to approximate it from the inverse-link.}
  \item{typical}{a function to be applied to the columns of the model matrix
    over which the effect is "averaged"; the default is \code{mean}.}
  \item{latent}{if \code{TRUE}, effects in a proportional-odds logit model
  	are computed on the scale of the latent response; if \code{FALSE} 
  	(the default) effects are computed as individual-level probabilities
  	and logits.}
  \item{x}{an object of class \code{"eff"} or \code{"effpoly"}.}
  \item{row.names, optional}{not used.}
  \item{...}{arguments to be passed down.}
}

\details{
  Normally, the functions to be used directly are \code{allEffects}, to return
  a list of high-order effects, and the generic \code{plot} function to plot the effects.
  (see \code{\link{plot.efflist}}, \code{\link{plot.eff}}, and \code{\link{plot.effpoly}}).
  Plots are drawn using the \code{\link[lattice]{xyplot}} (or in some cases, 
  the \code{\link[lattice]{densityplot}}) function in the 
  \code{lattice} package. Effects may also be printed (implicitly or explicitly via
  \code{print}) or summarized (using \code{summary})
  (see \code{\link{print.efflist}}, \code{\link{summary.efflist}},
  \code{\link{print.eff}}, \code{\link{summary.eff}}, \code{\link{print.effpoly}}, and \code{\link{summary.effpoly}}).
  
  If asked, the \code{effect} function will compute effects for terms that have 
  higher-order relatives in the model, averaging over those terms (which rarely makes sense), or for terms that
  do not appear in the model but are higher-order relatives of terms that do. 
  For example, for the model \code{Y ~ A*B + A*C + B*C}, one could
  compute the effect corresponding to the absent term \code{A:B:C}, which absorbs the constant, the
  \code{A}, \code{B}, and \code{C} main effects, and the three two-way interactions. In either of these
  cases, a warning is printed.
  
  In calculating effects, the strategy for `safe' prediction described
  in Hastie (1992: Sec. 7.3.3) is employed.
}

\value{
  For \code{lm} and \code{glm}, \code{effect} returns an \code{"eff"} object, and for \code{multinom}
  and \code{polr}, an \code{"effpoly"} object, with the following components:
  \item{term}{the term to which the effect pertains.}
  \item{formula}{the complete model formula.}
  \item{response}{a character string giving the name of the response variable.}
  \item{y.levels}{(for \code{"effpoly"} objects) levels of the polytomous response variable.}
  \item{variables}{a list with information about each predictor, including
    its name, whether it is a factor, and its levels or values.}
  \item{fit}{(for \code{"eff"} objects) a one-column matrix of fitted values, representing the effect
    on the scale of the linear predictor; this is a ravelled table, representing all
    combinations of predictor values.}
  \item{prob}{(for \code{"effpoly"} objects) a matrix giving fitted probabilities for the effect
  	for the various levels of the the response (columns) and combinations of the focal predictors (rows).}
  \item{logit}{(for \code{"effpoly"} objects) a matrix giving fitted logits for the effect
  	for the various levels of the the response (columns) and combinations of the focal predictors (rows).}
  \item{x}{a data frame, the columns of which are the predictors in the effect, and the
    rows of which give all combinations of values of these predictors.}
  \item{model.matrix}{the model matrix from which the effect was calculated.}
  \item{data}{a data frame with the data on which the fitted model was based.}
  \item{discrepancy}{the percentage discrepancy for the `safe' predictions of the original fit; should be very close to 0.}
  \item{se}{(for \code{"eff"} objects) a vector of standard errors for the effect, on the scale of the linear predictor.}
  \item{se.prob, se.logit}{(for \code{"effpoly"} objects) matrices of standard errors for the effect, on the probability and logit scales.}
  \item{lower, upper}{(for \code{"eff"} objects) one-column matrices of confidence limits, on the
    scale of the linear predictor.}
  \item{lower.prob, upper.prob, lower.logit, upper.logit}{(for \code{"effpoly"} objects) matrices of confidence limits
  	for the fitted logits and probabilities; the latter are computed by transforming
  	the former.}
  \item{confidence.level}{for the confidence limits.}
  \item{transformation}{(for \code{"eff"} objects) a two-element list, with element \code{link} giving the
    link function, and element \code{inverse} giving the inverse-link (mean) function.}
    
  \code{effectList} returns a list of \code{"eff"} or \code{"effpoly"} objects
  corresponding to the high-order terms of the model.
}

\references{
    Fox, J. (1987)
    Effect displays for generalized linear models.
    \emph{Sociological Methodology}
    \bold{17}, 347--361.
    
    Fox, J. (2003)
    Effect displays in R for generalised linear models.
    \emph{Journal of Statistical Software}
    \bold{8:15}, 1--27, <\url{http://www.jstatsoft.org/counter.php?id=75&url=v08/i15/effect-displays-revised.pdf&ct=1}>.
    
  	Fox, J. and R. Andersen (2006)
   	Effect displays  for multinomial and proportional-odds logit models.
  	\emph{Sociological Methodology} 
  	\bold{36}, 225--255.
    
    Hastie, T. J. (1992)
    Generalized additive models.
    In Chambers, J. M., and Hastie, T. J. (eds.)
    \emph{Statistical Models in S}, Wadsworth.
    }
    
    
    
\author{John Fox \email{jfox@mcmaster.ca} and Jangman Hong.}


\seealso{\code{\link{print.eff}}, \code{\link{summary.eff}}, \code{\link{plot.eff}}, 
    \code{\link{print.summary.eff}}, 
    \code{\link{print.effpoly}}, \code{\link{summary.effpoly}}, \code{\link{plot.effpoly}},
    \code{\link{print.efflist}}, \code{\link{summary.efflist}},
    \code{\link{plot.efflist}}, \code{\link[lattice]{xyplot}}, 
    \code{\link[lattice]{densityplot}}}

\examples{
mod.cowles <- glm(volunteer ~ sex + neuroticism*extraversion, 
    data=Cowles, family=binomial)
eff.cowles <- allEffects(mod.cowles, xlevels=list(neuroticism=0:24, 
    extraversion=seq(0, 24, 6)), given.values=c(sexmale=0.5))
eff.cowles

plot(eff.cowles, 'sex', ylab="Prob(Volunteer)")

plot(eff.cowles, 'neuroticism:extraversion', ylab="Prob(Volunteer)",
    ticks=list(at=c(.1,.25,.5,.75,.9)))

plot(eff.cowles, 'neuroticism:extraversion', multiline=TRUE, 
    ylab="Prob(Volunteer)")
    
plot(effect('sex:neuroticism:extraversion', mod.cowles,
    xlevels=list(neuroticism=0:24, extraversion=seq(0, 24, 6))), multiline=TRUE)

mod.beps <- multinom(vote ~ age + gender + economic.cond.national +
				economic.cond.household + Blair + Hague + Kennedy +
				Europe*political.knowledge, data=BEPS)
plot(effect("Europe*political.knowledge", mod.beps, 
	xlevels=list(Europe=1:11, political.knowledge=0:3)))
	
plot(effect("Europe*political.knowledge", mod.beps, 
		xlevels=list(Europe=1:11, political.knowledge=0:3),
		given.values=c(gendermale=0.5)),
	style="stacked", colors=c("blue", "red", "orange"), rug=FALSE)
	
mod.wvs <- polr(poverty ~ gender + religion + degree + country*poly(age,3),
	data=WVS)
plot(effect("country*poly(age, 3)", mod.wvs))

plot(effect("country*poly(age, 3)", mod.wvs), style="stacked", 
	colors=c("gray75", "gray50", "gray25"))
	
plot(effect("country*poly(age, 3)", latent=TRUE, mod.wvs))

mod.pres <- lm(prestige ~ log(income, 10) + poly(education, 3) + poly(women, 2), 
    data=Prestige)
eff.pres <- allEffects(mod.pres, default.levels=50)
plot(eff.pres, ask=FALSE)

}

\keyword{hplot}
\keyword{models}
