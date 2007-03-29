\name{effect}
\alias{effect}
\alias{effect.lm}
\alias{all.effects}
\alias{as.data.frame.eff}

\title{Functions For Constructing Effect Plots}
\description{
  \code{effect} constructs an \code{"effect"} object for a term (usually a high-order term) 
  in a linear or generalized linear model, absorbing the lower-order terms marginal
  to the term in question, and averaging over other terms in the model. 
  
  \code{all.effects} identifies all of the high-order terms in a model and returns
  a list of \code{"effect"} objects (i.e., an object of type \code{"effect.list"}).
}

\usage{
effect(term, mod, ...)

\method{effect}{lm}(term, mod, xlevels=list(), default.levels=10, se=TRUE, 
    confidence.level=.95, 
    transformation=list(link=family(mod)$linkfun, inverse=family(mod)$linkinv), 
    typical=mean, ...)
    
all.effects(mod, ...)

\method{as.data.frame}{eff}(x, row.names=NULL, optional=TRUE, ...)
}

\arguments{
  \item{term}{the quoted name of a term, usually, but not necessarily, a high-order 
    term in the model.}
  \item{mod}{an object of class \code{"lm"} or \code{"glm"}.}
  \item{xlevels}{an optional list of values at which to set covariates,
    with components of the form \code{covariate.name = vector.of.values}.}
  \item{default.levels}{number of values for covariates that are not
    specified explicitly via \code{xlevels}; covariate values set by
    default are evenly spaced between the minimum and maximum values in
    the data.}
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
  \item{...}{arguments to be passed down.}
  \item{x}{an object of type \code{"effect"}.}
  \item{row.names, optional}{not used.}
}

\details{
  Normally, the functions to be used directly are \code{all.effects}, to return
  a list of high-order effects, and the generic \code{plot} function to plot the effects.
  (see \code{\link{plot.eff.list}} and \code{\link{plot.eff}}).
  Plots are drawn using the \code{\link[lattice]{xyplot}} function in the 
  \code{lattice} package. Effects may also be printed (implicitly or explicitly via
  \code{print}) or summarized (using \code{summary})
  (see \code{\link{print.eff.list}}, \code{\link{summary.eff.list}},
  \code{\link{print.eff}}, and \code{\link{summary.eff}}).
  
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
  \code{effect} returns an \code{"eff"} object with the following components:
  \item{term}{the term to which the effect pertains.}
  \item{formula}{the complete model formula.}
  \item{response}{a character string giving the response variable.}
  \item{variables}{a list with information about each predictor, including
    its name, whether it is a factor, and its levels or values.}
  \item{fit}{a one-column matrix of fitted values, representing the effect
    on the scale of the linear predictor; this is a ravelled table, representing all
    combinations of predictor values.}
  \item{x}{a data frame, the columns of which are the predictors, and the
    rows of which give all combinations of values of the predictors.}
  \item{model.matrix}{the model matrix from which the effect was calculated.}
  \item{data}{a data frame with the data on which the fitted model was based.}
  \item{discrepancy}{the percentage discrepancy for the `safe' predictions of
    the original fit; should be very close to 0.}
  \item{se}{a vector of standard errors for the effect, on the scale of the linear predictor.}
  \item{lower, upper}{one-column matrices of confidence limits, on the
    scale of the linear predictor.}
  \item{confidence.level}{corresponding to the confidence limits.}
  \item{transformation}{a two-element list, with element \code{link} giving the
    link function, and element \code{inverse} giving the inverse-link (mean) function.}
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
    
    Hastie, T. J. (1992)
    Generalized additive models.
    In Chambers, J. M., and Hastie, T. J. (eds.)
    \emph{Statistical Models in S}, Wadsworth.
    }
    
    
    
\author{John Fox \email{jfox@mcmaster.ca}.}


\seealso{\code{\link{print.eff}}, \code{\link{summary.eff}}, \code{\link{plot.eff}}, 
    \code{\link{print.summary.eff}}, \code{\link{print.eff.list}}, \code{\link{summary.eff.list}},
    \code{\link{plot.eff.list}}, \code{\link[lattice]{xyplot}}}

\examples{
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

mod.pres <- lm(prestige ~ log(income, 10) + poly(education, 3) + poly(women, 2), 
    data=Prestige)
eff.pres <- all.effects(mod.pres, default.levels=50)
plot(eff.pres, ask=FALSE)

}

\keyword{hplot}
\keyword{models}
