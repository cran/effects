\name{effect}
\alias{effect}
\alias{Effect}
\alias{Effect.lm}
\alias{effect.lm}
\alias{effect.gls}
\alias{effect.multinom}
\alias{effect.polr}
\alias{effect.mer}
\alias{effect.lme}
\alias{allEffects}
\alias{allEffects.mer}
\alias{allEffects.lme}
\alias{as.data.frame.eff}
\alias{as.data.frame.effpoly}
\alias{as.data.frame.efflatent}
\alias{vcov.eff}

\title{Functions For Constructing Effect Plots}
\description{
  \code{effect} constructs an \code{"eff"} object for a term (usually a high-order term) 
  in a linear (fit by \code{lm} or \code{gls}) or generalized linear model (fit by \code{glm}), or an \code{"effpoly"} object for a term in a
   multinomial or proportional-odds logit model (fit respectively by \code{multinom} or \code{polr}), 
  absorbing the lower-order terms marginal
  to the term in question, and averaging over other terms in the model.  The function can
  also be used with mixed-effects models fit with \code{lmer} from the \code{lme4} package,
  or fit with \code{lme} from the \code{nlme} package.  In mixed effects models the
  analysis is for the fixed effects only, not for random effects. 
  
  \code{Effect} also constructs an an \code{"eff"} object, but rather than focusing on a term in the model,
  it focuses on a subset of the predictors. \code{Effect} is consequently more flexible and robust than 
  \code{effect}, and will work with some models for which \code{effect} fails, such as models with
  nested terms (see the examples).
  
  \code{allEffects} identifies all of the high-order terms in a model and returns
  a list of \code{"eff"} or \code{"effpoly"} objects (i.e., an object of type \code{"efflist"}).
}

\usage{
effect(term, mod, ...)

\method{effect}{lm}(term, mod, xlevels=list(), default.levels=10, given.values, 
	se=TRUE, confidence.level=.95, 
    transformation=list(link=family(mod)$linkfun, inverse=family(mod)$linkinv), 
    typical=mean, offset=mean, ...)
    
\method{effect}{gls}(term, mod, xlevels=list(), default.levels=10, given.values,
	se=TRUE, confidence.level=.95, transformation=NULL, typical=mean, ...)
    
\method{effect}{multinom}(term, mod, confidence.level=.95, xlevels=list(), default.levels=10, 
	given.values, se=TRUE, typical=mean, ...)
		
\method{effect}{polr}(term, mod, confidence.level=.95, xlevels=list(), default.levels=10, 
	given.values, se=TRUE, typical=mean, latent=FALSE, ...)
	
\method{effect}{mer}(term, mod, ...)

\method{effect}{lme}(term, mod, ...)

Effect(focal.predictors, mod, ...)

\method{Effect}{lm}(focal.predictors, mod, xlevels=list(), default.levels=10, given.values, 
		se=TRUE, confidence.level=0.95, 
		transformation=list(link=family(mod)$linkfun, inverse=family(mod)$linkinv), 
		typical=mean, offset=mean, ...)
    
allEffects(mod, ...)

\method{allEffects}{mer}(mod, ...)

\method{allEffects}{lme}(mod, ...)

\method{as.data.frame}{eff}(x, row.names=NULL, optional=TRUE, ...)

\method{as.data.frame}{effpoly}(x, row.names=NULL, optional=TRUE, ...)

\method{as.data.frame}{efflatent}(x, row.names=NULL, optional=TRUE, ...)

\method{vcov}{eff}(object, ...)
}

\arguments{
  \item{term}{the quoted name of a term, usually, but not necessarily, a high-order 
    term in the model. The term must be given exactly as it appears in the printed
    model, although either colons (\code{:}) or asterisks (\code{*}) may be used
    for interactions.}
  \item{focal.predictors}{a character vector of one or more predictors in the model.}
  \item{mod}{an object of class \code{"lm"}, \code{"gls"}, \code{"glm"}, 
    \code{"multinom"}, \code{"polr"} \code{"mer"} or \code{"lme"}. }
  \item{xlevels}{an optional list of values at which to set covariates,
    with components of the form \code{covariate.name = vector.of.values}.}
  \item{default.levels}{number of values for covariates that are not
    specified explicitly via \code{xlevels}; covariate values set by
    default are evenly spaced between the minimum and maximum values in
    the data.}
  \item{given.values}{a numeric vector of named elements, setting particular
  	columns of the model matrix to specific values for terms \emph{not}
  	appearing in an effect; if specified, this argument takes precedence over the
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
  \item{offset}{a function to be applied to the offset values (if there is an
    offset) in a linear or generalized linear model, or a mixed-effects
    model fit by \code{lmer} or \code{glmer}; or a numeric value, to which
    the offset will be set. The default is the \code{mean} function, and thus
    the offset will be set to its mean. \emph{Note:} Only offsets defined
    by the \code{offset} argument to \code{lm}, \code{glm}, \code{lmer}, or \code{glmer} will be handled
    correctly; use of the \code{offset} function in the model formula is
    not supported.}
  \item{latent}{if \code{TRUE}, effects in a proportional-odds logit model
  	are computed on the scale of the latent response; if \code{FALSE} 
  	(the default) effects are computed as individual-level probabilities
  	and logits.}
  \item{x}{an object of class \code{"eff"} or \code{"effpoly"}.}
  \item{row.names, optional}{not used.}
  \item{object}{an object of class \code{"eff"} for which the covariance matrix of the effects is desired.}
  \item{...}{arguments to be passed down.}
}

\details{
  Normally, the functions to be used directly are \code{allEffects}, to return
  a list of high-order effects, and the generic \code{plot} function to plot the effects.
  (see \code{\link{plot.efflist}}, \code{\link{plot.eff}}, and \code{\link{plot.effpoly}}).
  Alternatively, \code{Effect} can be used to vary a subset of predictors over their ranges,
  while other predictors are held to typical values.
  Plots are drawn using the \code{\link{xyplot}} (or in some cases, 
  the \code{\link{densityplot}}) function in the 
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
}

\value{
  For \code{lm}, \code{glm}, \code{mer} and \code{lme}, \code{effect} and \code{Effect} return 
  an \code{"eff"} object, and for \code{multinom}
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
  \item{offset}{value to which the offset is fixed; \code{0} if there is no offset.}
  \item{model}{(for \code{"effpoly"} objects) \code{"multinom"} or \code{"polr"}, as appropriate.}
  \item{vcov}{(for \code{"eff"} objects) a covariance matrix for the effect, on the scale of the linear predictor.}
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

\section{Warnings and Limitations}{
  The \code{effect} function handles factors and covariates differently, and becomes confused if one is changed to the other
  in a model formula. Consequently, formulas that include calls to \code{as.factor}, \code{factor}, or \code{numeric}
  (as, e.g., in \code{as.factor(income)}) will cause errors. Instead, create the modified variables outside of the model
  formula (e.g., \code{fincome <- as.factor(income)}) and use these in the model formula.
  
  Factors cannot have colons in level names (e.g., \code{"level:A"}); the \code{effect} function will confuse the
  colons with interactions; rename levels to remove or replace the colons (e.g., \code{"level.A"}).
  
  Offsets in linear and generalized linear models are supported, as are offsets in mixed models fit by 
  \code{lmer} or \code{glmer}, but must be supplied through the \code{offset}
  argument to \code{lm}, \code{glm}, \code{lmer} or \code{glmer}; 
  offsets supplied via calls to the \code{offset} function on the right-hand side
  of the model formula are not supported.
  }

\references{
    Fox, J. (1987)
    Effect displays for generalized linear models.
    \emph{Sociological Methodology}
    \bold{17}, 347--361.
    
    Fox, J. (2003)
    Effect displays in R for generalised linear models.
    \emph{Journal of Statistical Software}
    \bold{8:15}, 1--27, <\url{http://www.jstatsoft.org/v08/i15/}>.
    
  	Fox, J. and R. Andersen (2006)
   	Effect displays for multinomial and proportional-odds logit models.
  	\emph{Sociological Methodology} 
  	\bold{36}, 225--255.
  	
  	Fox, J. and J. Hong (2009).
    Effect displays in R for multinomial and proportional-odds logit models: 
    	    Extensions to the effects package.
    \emph{Journal of Statistical Software}
    \bold{32:1}, 1--24, <\url{http://www.jstatsoft.org/v32/i01/}>.
    
    Hastie, T. J. (1992)
    Generalized additive models.
    In Chambers, J. M., and Hastie, T. J. (eds.)
    \emph{Statistical Models in S}, Wadsworth.
    }
    
    
    
\author{John Fox \email{jfox@mcmaster.ca} and Jangman Hong.  Extension to
\code{mer} and \code{lme} objects by Sanford Weisberg \email{sandy@umn.edu}.}


\seealso{\code{\link{print.eff}}, \code{\link{summary.eff}}, \code{\link{plot.eff}}, 
    \code{\link{print.summary.eff}}, 
    \code{\link{print.effpoly}}, \code{\link{summary.effpoly}}, \code{\link{plot.effpoly}},
    \code{\link{print.efflist}}, \code{\link{summary.efflist}},
    \code{\link{plot.efflist}}, \code{\link{xyplot}}, 
    \code{\link{densityplot}}}

\examples{
mod.cowles <- glm(volunteer ~ sex + neuroticism*extraversion, 
    data=Cowles, family=binomial)
eff.cowles <- allEffects(mod.cowles, xlevels=list(neuroticism=0:24, 
    extraversion=seq(0, 24, 6)), given.values=c(sexmale=0.5))
eff.cowles

	# the following are equivalent:

eff.ne <- effect("neuroticism*extraversion", mod.cowles)
Eff.ne <- Effect(c("neuroticism", "extraversion"), mod.cowles)
all.equal(eff.ne$fit, Eff.ne$fit)

plot(eff.cowles, 'sex', ylab="Prob(Volunteer)")

plot(eff.cowles, 'neuroticism:extraversion', ylab="Prob(Volunteer)",
    ticks=list(at=c(.1,.25,.5,.75,.9)))

plot(eff.cowles, 'neuroticism:extraversion', multiline=TRUE, 
    ylab="Prob(Volunteer)")
    
plot(effect('sex:neuroticism:extraversion', mod.cowles,
    xlevels=list(neuroticism=0:24, extraversion=seq(0, 24, 6))), multiline=TRUE)
    
# a nested model:

mod <- lm(log(prestige) ~ income:type + education, data=Prestige)

# does not work: effect("income:type", mod, transformation=list(link=log, inverse=exp)) 

plot(Effect(c("income", "type"), mod, transformation=list(link=log, inverse=exp), 
	ylab="prestige")) # works
    
\donttest{ 
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
plot(effect("country*poly(age, 3)", mod.wvs), style="stacked")
plot(effect("country*poly(age, 3)", latent=TRUE, mod.wvs))
}

mod.pres <- lm(prestige ~ log(income, 10) + poly(education, 3) + poly(women, 2), 
    data=Prestige)
eff.pres <- allEffects(mod.pres, default.levels=50)
plot(eff.pres)

\dontrun{
	library(nlme) # for gls()
	mod.hart <- gls(fconvict ~ mconvict + tfr + partic + degrees, data=Hartnagel, 
	    correlation=corARMA(p=2, q=0), method="ML")
	plot(allEffects(mod.hart))

	library(lme4)
	data(cake, package="lme4")
	fm1 <- lmer(angle ~ recipe * temperature + (1|recipe:replicate), cake,
                   REML = FALSE)
	plot(effect("recipe:temperature", fm1), grid=TRUE)

	detach(package:lme4) # if previously attached
	library(nlme)
	data(cake, package="lme4")
	cake$rep <- with(cake, paste( as.character(recipe), as.character(replicate), sep=""))
	fm2 <- lme(angle ~ recipe * temperature, data=cake,
               random = ~ 1 | rep, method="ML")
	plot(effect("recipe:temperature", fm2), grid=TRUE)
}
}

\keyword{hplot}
\keyword{models}
