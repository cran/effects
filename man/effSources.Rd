\name{EffectMethods}
\alias{effSources}
\alias{effSources.default}
\alias{effSources.gls}
\alias{effSources.glmmPQL}
\alias{effSources.rlmerMod}
\alias{effSources.clm}
\alias{effSources.clm2}
\alias{effSources.clmm}
\alias{effSources.betareg}


\title{Functions For Constructing Effect Displays for Many Modeling Paradigms}

\description{The \pkg{effects} package can be used with many different classes of regression models beyond the basic models described in the \code{\link{Effect}} help page.  The \code{effSources} function is used to tell an \code{Effect} method where to find the information to compute effects.  Details are given in the vignette \href{../doc/methods-supported-by-effects.pdf}{Regression Models Supported by the effects Package}.  Unless you are trying to extend \pkg{effects} to a different class of regression model, you will have no use for this function.  
}

\usage{

\method{effSources}{default}(mod)

\method{effSources}{gls}(mod)

\method{effSources}{glmmPQL}(mod)

\method{effSources}{rlmerMod}(mod)

\method{effSources}{clm}(mod)

\method{effSources}{clm2}(mod)

\method{effSources}{clmm}(mod)

\method{effSources}{betareg}(mod)

}

\arguments{
  \item{mod}{a fitted model object of the appropriate class. The default method works for many regression models for which there is no specific method.}
  }

\details{
\code{Effect} function methods by default expect regression objects to provide certain quantities, including a model formula, coefficient estimates, and the variance matrix of the estimates, in standard places. \code{effSources} methods return a list of the necessary quantities if they are not in standard places.}

\value{See the vignette \href{../doc/methods-supported-by-effects.pdf}{Regression Models Supported by the effects Package}.}




\author{John Fox \email{jfox@mcmaster.ca}, Sanford Weisberg \email{sandy@umn.edu}}


%\examples{
%}

\keyword{internal}
