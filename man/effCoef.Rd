\name{effCoef}
\alias{effCoef}
\alias{effCoef.default}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to get coefficient estimates from regression models for use in the effects package.
}
\description{
This function uses the \code{\link[insight]{get_parameters}} function in the \code{insight} package to get a vector of regression coefficients for use in the effects package.  It converts the two-column \code{data.frame} returned by \code{get_parameters} to a vector of named elements.
}
\usage{

effCoef(mod, ...)

\method{effCoef}{default}(mod, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mod}{
A model object with a linear predictor representing fixed effects.  
}
\item{...}{Additional parameter passed to \code{\link[insight]{get_parameters}}.}
}
\details{
The \code{\link[insight]{get_parameters}} function can be used to retrieve the coefficient estimates corresponding to a linear predictor for many regression models, and return them as a two column \code{data.frame}, with regressor names in the first column and estimates in the second column.  This function converts this output to a named vector as is expected by the \code{effects} package.
}
\value{
A vector of coefficient estimates
%\references{
%% ~put references to the literature/web site here ~
}
\author{Sanford Weisberg \email{sandy@umn.edu}}

\seealso{\code{\link{get_parameters}}, and vignette \href{../doc/methods-supported-by-effects.pdf}{Regression Models Supported by the effects Package}}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
m1 <- lm(prestige ~ type + income + education, Duncan)
effCoef(m1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{models}
