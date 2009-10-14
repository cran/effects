\name{effects-package}
\Rdversion{1.1}
\alias{effects-package}
\alias{effects}
\docType{package}
\title{
Effect Displays for Linear, Generalized Linear, Multinomial-Logit, and Proportional-Odds Logit Models
}
\description{

Graphical and tabular effect displays, e.g., of interactions, for linear
generalized linear, multinomial-logit, and proportional-odds logit models.

}
\details{
\tabular{ll}{
Package: \tab effects\cr
Version: \tab 2.0-10\cr
Date: \tab 2009/10/14\cr
Depends: \tab R (>= 2.4.0), lattice, grid, MASS, nnet, colorspace\cr
LazyLoad: \tab yes\cr
LazyData: \tab yes\cr
License: \tab GPL (>= 2)\cr
URL: \tab http://www.r-project.org, http://socserv.socsci.mcmaster.ca/jfox/\cr
}

Index:
\preformatted{
Arrests                 Arrests for Marijuana Possession
BEPS                    British Election Panel Study
Cowles                  Cowles and Davis's Data on Volunteering
Prestige                Prestige of Canadian Occupations
Titanic                 Survival of Passengers on the Titanic
WVS                     World Values Surveys
Wells                   Well Switching in Bangladesh
effect                  Functions For Constructing Effect Plots
effects-deprecated      Deprecated Function in effects Package
summary.eff             Summarizing, Printing, and Plotting Effects
}

This package creates effect displays for various kinds of models, as explained in the references.
Typical usage is \code{plot(allEffects(model))} or \code{plot(allEffects(model), ask=FALSE)},
where \code{model} is an appropriate fitted-model object.
Additional arguments to \code{allEffects} and \code{plot} can be used to customize the resulting
displays. The function \code{effect} can be employed to produce an effect display for a 
particular term in the model, or to which terms in the model are marginal.
See \code{?effect} and \code{?plot.eff} for details.
}

\author{
John Fox <jfox@mcmaster.ca> and Jangman Hong. We are grateful to Robert Andersen, David Firth, and
Michael Friendly, for various suggestions.

Maintainer: John Fox <jfox@mcmaster.ca>
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
}


\keyword{ package }

