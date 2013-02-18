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
Version: \tab 2.2-4\cr
Date: \tab 2013/01/30\cr
Depends: \tab R (>= 2.4.0), lattice, grid, MASS, nnet, colorspace\cr
LazyLoad: \tab yes\cr
LazyData: \tab yes\cr
License: \tab GPL (>= 2)\cr
URL: \tab http://www.r-project.org, http://socserv.socsci.mcmaster.ca/jfox/\cr
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
John Fox <jfox@mcmaster.ca>, Sanford Weisberg, and Jangman Hong. We are grateful to Robert Andersen, David Firth, and
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
    
    Fox, J. and J. Hong (2009).
    Effect displays in R for multinomial and proportional-odds logit models: 
    	    Extensions to the effects package.
    \emph{Journal of Statistical Software}
    \bold{32:1}, 1--24, <\url{http://www.jstatsoft.org/v32/i01/}>.
}


\keyword{ package }

