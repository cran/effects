\name{effect}
\alias{effect}
\alias{Effect}
\alias{Effect.default}
\alias{Effect.lm}
\alias{Effect.mer}
\alias{Effect.merMod}
\alias{Effect.lme}
\alias{Effect.gls}
\alias{Effect.multinom}
\alias{Effect.polr}
\alias{Effect.poLCA}
\alias{Effect.clm2}
\alias{Effect.clm}
\alias{Effect.clmm}
\alias{Effect.mlm}
\alias{Effect.svyglm}
\alias{allEffects}
\alias{as.data.frame.eff}
\alias{as.data.frame.effpoly}
\alias{as.data.frame.efflatent}
\alias{vcov.eff}

\title{Functions For Constructing Effect Displays}
\description{
  \code{Effect} and \code{effect} construct an \code{"eff"} object for a term (usually a high-order term)
  in a linear model (fit by \code{\link{lm}} or \code{\link[nlme]{gls}}) or generalized linear model (fit by \code{\link{glm}}), or an
  \code{"effpoly"} object for a term in a
  multinomial or proportional-odds logit model (fit respectively by \code{\link[nnet]{multinom}} or \code{\link[MASS]{polr}}),
  absorbing the lower-order terms marginal
  to the term in question, and averaging over other terms in the model.  For multivariate linear models
  (of class \code{"mlm"}, fit by \code{\link{lm}}), the function constructs a list of \code{"eff"} objects separately for the various
  response variables.

  \code{effect} builds the required object by specifying explicity a focal term like \code{"a:b"}
  for an \code{a} by \code{b} interaction.  \code{Effect} specifies the predictors in
  the term, for example \code{c("a", "b")}, rather than the term itself.  \code{Effect} is consequently
  more flexible and robust than
  \code{effect}, and will succeed with some models for which \code{effect} fails.  The \code{effect} function
  works by constructing a call to \code{Effect}.

  The \code{Effect} and \code{effect} functions can also be used with some mixed-effects models,
  including linear and generalized linear mixed-effects models
  fit by \code{\link[lme4]{lmer}} and \code{\link[lme4]{glmer}} from the \pkg{lme4} package
  and \code{\link[nlme]{lme}} from the \pkg{nlme} package; ordinal logit mixed effects fit with \code{\link[ordinal]{clm2}} or \code{\link[ordinal]{clmm}} from the
  \pkg{ordinal} package, polytomous latent-class models fit by the \code{\link[poLCA]{poLCA}} function in
  the \pkg{poLCA} package, and generalized linear models fit to data from complex surveys
  using the \code{\link[survey]{svyglm}} function in the \pkg{survey} package.  The displays for mixed-effects models are for the fixed effects only, not for random effects.

  There is a default
  method for \code{Effect} that should work with most model objects that have a linear predictor and that respond
  to the \code{\link{coef}}, \code{\link{model.frame}}, \code{\link{formula}}, and \code{\link{vcov}} functions.

  \code{allEffects} identifies all of the high-order terms in a model and returns
  a list of \code{"eff"} or \code{"effpoly"} objects (i.e., an object of type \code{"efflist"}).

  For information on computing and displaying \emph{predictor effects}, see \code{\link{predictorEffect}} and \code{\link{plot.predictoreff}}.

  For further information about plotting effects, see \code{\link{plot.eff}}.
}

\usage{
effect(term, mod, vcov.=vcov, ...)

Effect(focal.predictors, mod, ...)

\method{Effect}{lm}(focal.predictors, mod, xlevels=list(),  
    fixed.predictors, vcov. = vcov, confint=TRUE,
    transformation = list(link = family(mod)$linkfun, 
        inverse = family(mod)$linkinv), 
    partial.residuals=FALSE, quantiles=seq(0.2, 0.8, by=0.2),
    x.var=NULL,  ...,
    #legacy arguments:
    given.values, typical, offset, se, confidence.level)

\method{Effect}{gls}(focal.predictors, mod, xlevels = list(), 
    fixed.predictors, vcov. = vcov, confint=TRUE, 
    transformation = NULL, ...,
    #legacy arguments:
    given.values, typical, se, confidence.level)

\method{Effect}{multinom}(focal.predictors, mod, 
    xlevels=list(), fixed.predictors,
    vcov. = vcov, confint=TRUE, ...,
    #legacy arguments:
    se, confidence.level, given.values, typical)

\method{Effect}{polr}(focal.predictors, mod, 
    xlevels=list(), fixed.predictors,
    vcov.=vcov, confint=TRUE, latent=FALSE, ...,
    #legacy arguments:
    se, confidence.level, given.values, typical)

\method{Effect}{clm2}(focal.predictors, mod, ...)

\method{Effect}{clmm}(focal.predictors, mod, ...)

\method{Effect}{clm}(focal.predictors, mod, ...)

\method{Effect}{mer}(focal.predictors, mod,
    KR=FALSE, ...)

\method{Effect}{merMod}(focal.predictors, mod,
    KR=FALSE, ...)

\method{Effect}{lme}(focal.predictors, mod, ...)

\method{Effect}{poLCA}(focal.predictors, mod, ...)

\method{Effect}{mlm}(focal.predictors, mod, response, ...)

\method{Effect}{svyglm}(focal.predictors, mod, fixed.predictors, ...)

\method{Effect}{default}(focal.predictors, mod, xlevels = list(), 
    fixed.predictors, vcov. = vcov, confint=TRUE, 
    transformation = list(link = I, inverse = I), ...,
    #legacy arguments:
    se, confidence.level, given.values, typical, offset)

allEffects(mod, ...)


\method{as.data.frame}{eff}(x, row.names=NULL, optional=TRUE,
    transform=x$transformation$inverse, ...)

\method{as.data.frame}{effpoly}(x, row.names=NULL, optional=TRUE, ...)

\method{as.data.frame}{efflatent}(x, row.names=NULL, optional=TRUE, ...)

\method{vcov}{eff}(object, ...)
}

\arguments{
  \item{term}{the quoted name of a term, usually, but not necessarily, a high-order
    term in the model. The term must be given exactly as it appears in the printed
    model, although either colons (\code{:}) or asterisks (\code{*}) may be used
    for interactions.}
  \item{focal.predictors}{a character vector of one or more predictors in the model in any order.}
  \item{mod}{an object of class \code{"lm"}, \code{"gls"}, \code{"glm"},
    \code{"multinom"}, \code{"polr"}, \code{"mer"} (or \code{"merMod"}), \code{"lme"},
    or \code{"poLCA"}. }
  \item{xlevels}{this argument is used to set the number of levels for any
    focal predictor that is not a factor.  If \code{xlevels=NULL}, then
    each numeric predictor is represented by five values equally spaced over
    its range and then rounded to 'nice' numbers.  If \code{xlevels=n}
    is an integer, then each numeric predictor is represented by \code{n} equally
    spaced values rounded to 'nice' numbers.  More generally, \code{xlevels} can be a named list of values
    at which to set each numeric predictor.  For example, \code{xlevels=list(x1=c(2, 4, 7),
    x2=5)} would use the values 2, 4 and 7 for the levels of \code{x1}, use 5
    equally spaced levels for the levels of \code{x2}, and use the default for
    any other numeric predictors. If partial residuals are computed, then the focal
    predictor that is to appear on the horizontal axis of an effect plot is
    evaluated at 100 equally spaced values along its full range, and, by default,
    other numeric predictors are evaluated at the quantiles specified in the
    \code{quantiles} argument, unless their values are given explicitly in
    \code{xlevels}.}
  \item{fixed.predictors}{an optional list of specifications affecting the values at which fixed predictors
    for an effect are set, potentially including:
    
    \describe{
    
      \item{given.values}{a numeric vector of named elements, setting particular
  	columns of the model matrix to specific values for predictors that are \emph{not}
  	focal predictors; if specified, these takes precedence over the
  	application of the function given in the \code{typical} list element
  	(below). Care must be taken in specifying these values --- e.g.,
  	for a factor, the values of all contrasts should be given and these
  	should be consistent with each other.}
  	
  	  \item{typical}{a function to be applied to the columns of the model matrix
    over which the effect is "averaged"; with the exception of the \code{"svyglm"} method,
    the default is \code{\link{mean}}. For\code{"svyglm"} objects, the default is to
    use the survey-design weighted mean.}
    
      \item{apply.typical.to.factors}{It generally doesn't make sense to apply typical
    values that aren't means (e.g., medians) to the columns of the model-matrix
    representing contrasts for factors. This value generally defaults to \code{FALSE}
    except for \code{"svyglm"} objects, for which the default is \code{TRUE}, using the the survey-design
    weighted mean.}
  
      \item{offset}{a function to be applied to the offset values (if there is an
    offset) in a linear or generalized linear model, or a mixed-effects
    model fit by \code{\link[lme4]{lmer}} or \code{\link[lme4]{glmer}}; or a numeric value, to which
    the offset will be set. The default is the \code{\link{mean}} function, and thus
    the offset will be set to its mean; in the case of \code{"svyglm"}
    objects, the default is to use the survey-design weighted mean. 
    \emph{Note:} Only offsets defined
    by the \code{offset} argument to \code{\link{lm}}, \code{\link{glm}}, \code{\link[survey]{svyglm}}, 
    \code{\link[lme4]{lmer}}, or \code{\link[lme4]{glmer}} will be handled
    correctly; use of the \code{offset} function in the model formula is
    not supported.}
    
    }
  }

  \item{vcov.}{A function or the name of a function that will be used to get the estimated variance-covariance
    matrix of the estimated coefficients.  This will ordinarily be the default,
    \code{\link{vcov}}, which will result in the function call \code{vcov(mod)} to get the variance-covariance matrix.
    You can use the name of any function that takes the model object as
    its first argument and returns an estimated sample covariance matrix, such as the \code{\link[car]{hccm}}
    function in the \pkg{car} package, which returns a heteroscedasticity corrected estimate for a
    linear model.}
  \item{confint}{a list with any or all of the following elements, controlling whether and how standard errors
    and confidence limits are computed for the effects: \code{compute} (default \code{TRUE}), whether or not
    to compute standard errors and confidence limits; \code{level} (default \code{0.95}), confidence level for
    confidence limits; \code{type}, one of \code{"pointwise"} (the default), \code{"Scheffe"}, or \code{"scheffe"},
    whether to compute confidence limits with specified coverage at each point for an effect or to compute limits
    for a Scheffe-type confidence envelope. For \code{mer}, \code{merMod}, and \code{lme} objects, the
    normal distribution is used to get confidence limits.}
  \item{KR}{if \code{TRUE} and the \pkg{pbkrtest} package is installed,
    use the Kenward-Roger coefficient covariance
    matrix to compute effect standard errors for linear mixed models fit with \code{\link[lme4]{lmer}}
    in the \pkg{lme4} package.  The default is \code{FALSE} because the computation can be
    very slow.}
  \item{transformation}{a two-element list with elements \code{link} and
    \code{inverse}. For a generalized linear model, these are by default
    the link function and inverse-link (mean) function. For a linear model,
    these default to \code{NULL}. If \code{NULL}, the identity (or inhibit) function,
    \code{\link{I}}, is used; this effect can also be achieved by setting the
    argument to \code{NULL}. The inverse-link may be used to transform effects
    when they are printed or plotted; the link may be used in positioning
    axis labels (see below). If the link is not given, an attempt will be
    made to approximate it from the inverse-link.}
  \item{partial.residuals}{if \code{TRUE}, residuals for a linear or generalized linear model
    will be computed and saved; if \code{FALSE}
    (the default), residuals are suppressed. If residuals are saved, partial residuals are computed when
    the effect is plotted: see \code{\link{plot.eff}}. This argument may also be used for mixed-effects models.}
  \item{quantiles}{quantiles at which to evaluate numeric focal predictors \emph{not} on the
    horizontal axis, used only when partial residuals are displayed; superceded if the \code{xlevels} argument
    gives specific values for a predictor.}
  \item{x.var}{the name or index of the numeric predictor to define the horizontal axis of an effect
    plot for a linear or generalized linear model; the default is \code{NULL}, in which case the first
    numeric predictor in the effect will be used \emph{if} partial residuals are to be computed. This argument
    is intended to be used when \code{partial.residuals} is \code{TRUE}; otherwise, the variable on the horizontal
    axis can be chosen when the effect object is plotted: see \code{\link{plot.eff}}.}
  \item{latent}{if \code{TRUE}, effects in a proportional-odds logit model
  	are computed on the scale of the latent response; if \code{FALSE}
  	(the default) effects are computed as individual-level probabilities
  	and logits.}
  \item{x}{an object of class \code{"eff"}, \code{"effpoly"}, or \code{"efflatent"}.}
  \item{transform}{a transformation to be applied to the effects and confidence limits,
    by default taken from the inverse link function saved in the \code{"eff"} object.}
  \item{row.names, optional}{not used.}
  \item{response}{for an \code{"mlm"} object, a vector containing the name(s) or indices of one or more response variable(s).
    The default is to use all responses in the model.}
  \item{object}{an object of class \code{"eff"} for which the covariance matrix of the effects is desired.}
  \item{...}{arguments to be passed down.}
  \item{se, confidence.level, given.values, typical, offset}{legacy arguments retained for backwards compatability; if present, these arguments take precedence
  of the \code{compute} and \code{level} elements over the \code{confint} list argument and the
  \code{given.values}, \code{typical}, and \code{offset} elements of the \code{fixed.predictors}
  list argument. See \code{\link{LegacyArguments}} for details.}
}

\details{
  Normally, the functions to be used directly are \code{allEffects}, to return
  a list of high-order effects, and the generic \code{plot} function to plot the effects.
  (see \code{\link{plot.efflist}}, \code{\link{plot.eff}}, and \code{\link{plot.effpoly}}).
  Alternatively, \code{Effect} can be used to vary a subset of predictors over their ranges,
  while other predictors are held to typical values.
  Plots are drawn using the \code{\link{xyplot}} (or in some cases,
  the \code{\link{densityplot}}) function in the
  \pkg{lattice} package. Effects may also be printed (implicitly or explicitly via
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

  The \code{as.data.frame} methods convert effect objects to data frames to facilitate the construction
  of custom displays. In the case of \code{"eff"} objects, the \code{se} element in the data frame is always
  on the scale of the linear predictor, and the transformation used for the fit and confidence limits is saved in
  a \code{"transformation"} attribute.

  See \code{\link{predictorEffects}} for an alternative paradigm for getting effects.
}

\value{
  For \code{lm}, \code{glm},\code{svyglm}, \code{mer} and \code{lme}, \code{effect} and \code{Effect} return
  an \code{"eff"} object, and for \code{multinom},
  \code{polr}, \code{clm}, \code{clmm} and \code{clm2}, an \code{"effpoly"} object, with the components listed below.
  For an \code{"mlm"} object with one response specified, an \code{"eff"} object is returned, otherwise an \code{"efflist"} object
  is returned, containing one \code{"eff"} object for each \code{response}.
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
  \item{discrepancy}{the percentage discrepancy for the `safe' predictions of the original fit; should be very close to 0.
    Note: except for \code{gls} models, this is now necessarily 0.}
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
  \item{residuals}{(working) residuals for linear or generalized linear models, to be used by
  \code{\link{plot.eff}} to plot partial residuals.}
  \item{x.var}{the name of the predictor to appear on the horizontal axis of an effect plot made from the
    returned object; will usually be \code{NULL} if partial residuals aren't computed.}
  \item{family}{for a \code{"glm"} model, the name of the distributional family of the model;
    for an \code{"lm"} model, this is \code{"gaussian"}; otherwise \code{NULL}. The \code{family} controls
    how partial residuals are smoothed in plots.}

  \code{allEffects} returns an \code{"efflist"} object, a list of \code{"eff"} or \code{"effpoly"} objects
  corresponding to the high-order terms of the model.

If \code{mod} is of class \code{"poLCA"} (from the \code{poLCA} package), representing a
polytomous latent class model, effects are computed for the predictors given the
estimated latent classes.
The result is of class \code{"eff"} if the latent class model has 2 categories
and of class \code{"effpoly"} with more than 2 categories.
}

\section{Warnings and Limitations}{
  The \code{Effect} function handles factors and covariates differently, and is likely to be confused if one is changed to the other
  in a model formula. Consequently, formulas that include calls to \code{as.factor}, \code{factor}, or \code{numeric}
  (as, e.g., in \code{y ~ as.factor(income)}) will cause errors. Instead, create the modified variables outside of the model
  formula (e.g., \code{fincome <- as.factor(income)}) and use these in the model formula.

  Factors cannot have colons in level names (e.g., \code{"level:A"}); the \code{effect} function will confuse the
  colons with interactions; rename levels to remove or replace the colons (e.g., \code{"level.A"}).

  The functions in the \pkg{effects} package work properly with predictors that are numeric or factors;
  consequently, e.g., convert logical predictors to factors, and dates to numeric.

  Empty cells in crossed-factors are now permitted for \code{"lm"}, \code{"glm"}, and \code{"multinom"} models.
  For \code{"multinom"} models with two or more crossed factors with an empty cell, stacked area plots
  apparently do not work because of a bug in the \code{\link[lattice]{barchart}} function in the \pkg{lattice} package.  However, the default
  line plots do work.

  Offsets in linear and generalized linear models are supported, as are offsets in mixed models fit by
  \code{lmer} or \code{glmer}, but must be supplied through the \code{offset}
  argument to \code{lm}, \code{glm}, \code{lmer} or \code{glmer};
  offsets supplied via calls to the \code{offset} function on the right-hand side
  of the model formula are not supported.

  Fitting ordinal mixed-models using \code{\link[ordinal]{clmm}} or \code{\link[ordinal]{clmm2}} permits many options, including a variety of link functions,
  scale functions, nominal regressors, and various methods for setting thresholds.  Effects are currently generated
  only for the default values of the arguments \code{scale}, \code{nominal}, \code{link} and \code{threshold},
  which is equivalent to fitting an ordinal response mixed effects model with a logit link.
  The effect methods can also be
  used with objects created using \code{\link[ordinal]{clm}} or \code{\link[ordinal]{clm2}} fitting ordinal response logit models with no random effects, with
  results similar to those from \code{\link[MASS]{polr}} in the \pkg{MASS} package.

  Calling any of these functions from within a user-written function may result in errors due
  to R's scoping rules.  See the vignette \code{embedding.pdf} for the \pkg{car} package
  for a solution to this problem.
  }

\references{
    Fox, J. (1987).
    Effect displays for generalized linear models.
    \emph{Sociological Methodology}
    \bold{17}, 347--361.

    Fox, J. (2003)
    Effect displays in R for generalised linear models.
    \emph{Journal of Statistical Software}
    \bold{8:15}, 1--27, <\url{http://www.jstatsoft.org/v08/i15/}>.

  	Fox, J. and R. Andersen (2006).
   	Effect displays for multinomial and proportional-odds logit models.
  	\emph{Sociological Methodology}
  	\bold{36}, 225--255.

  	Fox, J. and J. Hong (2009).
    Effect displays in R for multinomial and proportional-odds logit models:
    	    Extensions to the effects package.
    \emph{Journal of Statistical Software}
    \bold{32:1}, 1--24, <\url{http://www.jstatsoft.org/v32/i01/}>.

    Fox, J. and S. Weisberg (forthcoming).
    Visualizing Fit and Lack of Fit in Complex Regression Models:
        Effect Plots with Partial Residuals.
    \emph{Journal of Statistical Software}.

    Hastie, T. J. (1992).
    Generalized additive models.
    In Chambers, J. M., and Hastie, T. J. (eds.)
    \emph{Statistical Models in S}, Wadsworth.

    Weisberg, S. (2014).
    \emph{Applied Linear Regression}, 4th edition, Wiley,
    \url{http://z.umn.edu/alr4ed}.
    }



\author{John Fox \email{jfox@mcmaster.ca}, Sanford Weisberg \email{sandy@umn.edu}
and Jangman Hong.}


\seealso{\code{\link{LegacyArguments}}. For information on printing, summarizing, and plotting effects:
    \code{\link{print.eff}}, \code{\link{summary.eff}}, \code{\link{plot.eff}},
    \code{\link{print.summary.eff}},
    \code{\link{print.effpoly}}, \code{\link{summary.effpoly}}, \code{\link{plot.effpoly}},
    \code{\link{print.efflist}}, \code{\link{summary.efflist}},
    \code{\link{plot.efflist}}, \code{\link{xyplot}},
    \code{\link{densityplot}}.}

\examples{
mod.cowles <- glm(volunteer ~ sex + neuroticism*extraversion,
                  data=Cowles, family=binomial)
eff.cowles <- allEffects(mod.cowles, xlevels=list(extraversion=seq(0, 24, 6)),
                   fixed.predictors=list(given.values=c(sexmale=0.5)))
eff.cowles
as.data.frame(eff.cowles[[2]])

\donttest{
# the following are equivalent:
eff.ne <- effect("neuroticism*extraversion", mod.cowles)
Eff.ne <- Effect(c("neuroticism", "extraversion"), mod.cowles)
all.equal(eff.ne$fit, Eff.ne$fit)

plot(eff.cowles, 'sex', axes=list(y=list(lab="Prob(Volunteer)")))

plot(eff.cowles, 'neuroticism:extraversion',
     axes=list(y=list(lab="Prob(Volunteer)",
        ticks=list(at=c(.1,.25,.5,.75,.9)))))
        
plot(Effect(c("neuroticism", "extraversion"), mod.cowles,
            confint=list(type="Scheffe"),
            xlevels=list(extraversion=seq(0, 24, 6)), 
            fixed.predictors=list(given.values=c(sexmale=0.5))),
     axes=list(y=list(lab="Prob(Volunteer)",
        ticks=list(at=c(.1,.25,.5,.75,.9)))))

plot(eff.cowles, 'neuroticism:extraversion', lines=list(multiline=TRUE),
     axes=list(y=list(lab="Prob(Volunteer)")))

plot(effect('sex:neuroticism:extraversion', mod.cowles,
            xlevels=list(extraversion=seq(0, 24, 6))),
     lines=list(multiline=TRUE))
}

# a nested model:

mod <- lm(log(prestige) ~ income:type + education, data=Prestige)

plot(Effect(c("income", "type"), mod, transformation=list(link=log, inverse=exp)),
     axes=list(y=list(lab="prestige")))


if (require(nnet)){
    mod.beps <- multinom(vote ~ age + gender + economic.cond.national +
                             economic.cond.household + Blair + Hague + Kennedy +
                             Europe*political.knowledge, data=BEPS)
    \donttest{
    plot(effect("Europe*political.knowledge", mod.beps,
                xlevels=list(political.knowledge=0:3)))
    }

    plot(Effect(c("Europe", "political.knowledge"), mod.beps,
                xlevels=list(Europe=1:11, political.knowledge=0:3),
                fixed.predictors=list(given.values=c(gendermale=0.5))),
         lines=list(col=c("blue", "red", "orange")),
         axes=list(x=list(rug=FALSE), y=list(style="stacked")))

    \donttest{
    plot(effect("Europe*political.knowledge", mod.beps, # equivalent
                xlevels=list(Europe=1:11, political.knowledge=0:3),
                fixed.predictors=list(given.values=c(gendermale=0.5))),
         lines=list(col=c("blue", "red", "orange")),
         axes=list(x=list(rug=FALSE), y=list(style="stacked")))
    }
}

if (require(MASS)){
    mod.wvs <- polr(poverty ~ gender + religion + degree + country*poly(age,3),
                    data=WVS)
    \donttest{
    plot(effect("country*poly(age, 3)", mod.wvs))
    }

    plot(Effect(c("country", "age"), mod.wvs),
         axes=list(y=list(style="stacked")))

    \donttest{
    plot(effect("country*poly(age, 3)", mod.wvs),
         axes=list(y=list(style="stacked"))) # equivalent

    plot(effect("country*poly(age, 3)", latent=TRUE, mod.wvs))
    plot(effect("country*poly(age, 3)", latent=TRUE, mod.wvs, 
         confint=list(type="scheffe"))) # Scheffe-type confidence envelopes
    }
}


mod.pres <- lm(prestige ~ log(income, 10) + poly(education, 3) + poly(women, 2),
               data=Prestige)
eff.pres <- allEffects(mod.pres, xlevels=50)
plot(eff.pres)
plot(eff.pres[1],
     axes=list(x=list(income=list(
             transform=list(trans=log10, inverse=function(x) 10^x),
             ticks=list(at=c(1000, 2000, 5000, 10000, 20000))
    ))))
\donttest{
# linear model with log-response and log-predictor
# to illustrate transforming axes and setting tick labels
mod.pres1 <- lm(log(prestige) ~ log(income) + poly(education, 3) + poly(women, 2),
                data=Prestige)
# effect of the log-predictor
eff.log <- Effect("income", mod.pres1)
# effect of the log-predictor transformed to the arithmetic scale
eff.trans <- Effect("income", mod.pres1, transformation=list(link=log, inverse=exp))
#variations:
# y-axis:  scale is log, tick labels are log
# x-axis:  scale is arithmetic, tick labels are arithmetic
plot(eff.log)

# y-axis:  scale is log, tick labels are log
# x-axis:  scale is log, tick labels are arithmetic
plot(eff.log, axes=list(x=list(income=list(
    transform=list(trans=log, inverse=exp),
    ticks=list(at=c(5000, 10000, 20000)),
    lab="income, log-scale"))))

# y-axis:  scale is log, tick labels are arithmetic
# x-axis:  scale is arithmetic, tick labels are arithmetic
plot(eff.trans, axes=list(y=list(lab="prestige")))

# y-axis:  scale is arithmetic, tick labels are arithmetic
# x-axis:  scale is arithmetic, tick labels are arithmetic
plot(eff.trans, axes=list(y=list(type="response", lab="prestige")))

# y-axis:  scale is log, tick labels are arithmetic
# x-axis:  scale is log, tick labels are arithmetic
plot(eff.trans, axes=list(
       x=list(income=list(
            transform=list(trans=log, inverse=exp),
            ticks=list(at=c(1000, 2000, 5000, 10000, 20000)),
            lab="income, log-scale")),
        y=list(lab="prestige, log-scale")),
     main="Both response and X in log-scale")

# y-axis:  scale is arithmetic, tick labels are arithmetic
# x-axis:  scale is log, tick labels are arithmetic
plot(eff.trans, axes=list(
        x=list(
            income=list(transform=list(trans=log, inverse=exp),
                        ticks=list(at=c(1000, 2000, 5000, 10000, 20000)),
                        lab="income, log-scale")),
        y=list(type="response", lab="prestige")))
}

if (require(nlme)){ # for gls()
    mod.hart <- gls(fconvict ~ mconvict + tfr + partic + degrees, data=Hartnagel,
                    correlation=corARMA(p=2, q=0), method="ML")
    plot(allEffects(mod.hart))
    detach(package:nlme)
}

if (require(lme4)){
    data(cake, package="lme4")
    fm1 <- lmer(angle ~ recipe * temperature + (1|recipe:replicate), cake,
                REML = FALSE)
    plot(Effect(c("recipe", "temperature"), fm1))
    \donttest{
    plot(effect("recipe:temperature", fm1),
         axes=list(grid=TRUE)) # equivalent (plus grid)
    }
    if (any(grepl("pbkrtest", search()))) detach(package:pbkrtest)
    detach(package:lme4)
}

\donttest{
if (require(nlme) && length(find.package("lme4", quiet=TRUE)) > 0){
    data(cake, package="lme4")
    cake$rep <- with(cake, paste( as.character(recipe), as.character(replicate), sep=""))
    fm2 <- lme(angle ~ recipe * temperature, data=cake,
               random = ~ 1 | rep, method="ML")
    plot(Effect(c("recipe", "temperature"), fm2))
    plot(effect("recipe:temperature", fm2),
         axes=list(grid=TRUE))  # equivalent (plus grid)
    }
    detach(package:nlme)
}

\donttest{
if (require(poLCA)){
    data(election)
    f2a <- cbind(MORALG,CARESG,KNOWG,LEADG,DISHONG,INTELG,
                 MORALB,CARESB,KNOWB,LEADB,DISHONB,INTELB)~PARTY*AGE
    nes2a <- poLCA(f2a,election,nclass=3,nrep=5)
    plot(Effect(c("PARTY", "AGE"), nes2a),
         axes=list(y=list(style="stacked")))
}
}

# mlm example
if (require(heplots)) {
    data(NLSY, package="heplots")
    mod <- lm(cbind(read,math) ~ income+educ, data=NLSY)
    eff.inc <- Effect("income", mod)
    plot(eff.inc)
    eff.edu <- Effect("educ", mod)
    plot(eff.edu, axes=list(x=list(rug=FALSE), grid=TRUE))
    \donttest{
    plot(Effect("educ", mod, response="read"))
    }
    detach(package:heplots)
}

# svyglm() example (adapting an example from the survey package)

\donttest{
if (require(survey)){
  data(api)
  dstrat<-svydesign(id=~1, strata=~stype, weights=~pw, 
    data=apistrat, fpc=~fpc) 
  mod <- svyglm(sch.wide ~ ell + meals + mobility, design=dstrat, 
    family=quasibinomial())
  plot(allEffects(mod), 
    axes=list(y=list(lim=log(c(0.4, 0.99)/c(0.6, 0.01)), 
      ticks=list(at=c(0.4, 0.75, 0.9, 0.95, 0.99)))))
}
}

# component + residual plot examples

\donttest{
Prestige$type <- factor(Prestige$type, levels=c("bc", "wc", "prof"))

mod.prestige.1 <- lm(prestige ~ income + education, data=Prestige)
plot(allEffects(mod.prestige.1, partial.residuals=TRUE)) # standard C+R plots
plot(allEffects(mod.prestige.1, partial.residuals=TRUE,
    confint=list(type="scheffe"))) # with Scheffe-type confidence bands

mod.prestige.2 <- lm(prestige ~ type*(income + education), data=Prestige)
plot(allEffects(mod.prestige.2, partial.residuals=TRUE))

mod.prestige.3 <- lm(prestige ~ type + income*education, data=Prestige)
plot(Effect(c("income", "education"), mod.prestige.3, partial.residuals=TRUE),
     partial.residuals=list(span=1))
}

#  artificial data

set.seed(12345)
x1 <- runif(500, -75, 100)
x2 <- runif(500, -75, 100)
y <- 10 + 5*x1 + 5*x2 + x1^2 + x2^2 + x1*x2 + rnorm(500, 0, 1e3)
Data <- data.frame(y, x1, x2)
mod.1 <- lm(y ~ poly(x1, x2, degree=2, raw=TRUE), data=Data)
# raw=TRUE necessary for safe prediction
mod.2 <- lm(y ~ x1*x2, data=Data)
mod.3 <- lm(y ~ x1 + x2, data=Data)

plot(Effect(c("x1", "x2"), mod.1, partial.residuals=TRUE)) # correct model
plot(Effect(c("x1", "x2"), mod.2, partial.residuals=TRUE)) # wrong model
plot(Effect(c("x1", "x2"), mod.3, partial.residuals=TRUE)) # wrong model
}

\keyword{hplot}
\keyword{models}
