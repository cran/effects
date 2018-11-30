\name{summary.eff}
\alias{print.eff}
\alias{print.effpoly}
\alias{print.efflatent}
\alias{print.efflist}
\alias{print.mlm.efflist}
\alias{print.summary.eff}
\alias{summary.eff}
\alias{summary.effpoly}
\alias{summary.efflatent}
\alias{summary.efflist}
\alias{summary.mlm.efflist}
\alias{as.data.frame.eff}
\alias{as.data.frame.effpoly}
\alias{as.data.frame.efflatent}
\alias{as.data.frame.efflist}
\alias{vcov.eff}

\title{Summarizing and Printing Effects}

\description{
\code{summary}, \code{print}, and \code{as.data.frame} methods for objects created using the effects package.
}

\usage{
\method{print}{eff}(x, type=c("response", "link"), ...)
\method{print}{effpoly}(x, type=c("probability", "logits"), ...)
\method{print}{efflatent}(x, ...)
\method{print}{efflist}(x, ...)
\method{print}{mlm.efflist}(x, ...)
\method{print}{summary.eff}(x, ...)
\method{summary}{eff}(object, type=c("response", "link"), ...)
\method{summary}{effpoly}(object, type=c("probability", "logits"), ...)
\method{summary}{efflatent}(object, ...)
\method{summary}{efflist}(object, ...)
\method{summary}{mlm.efflist}(object, ...)
\method{as.data.frame}{eff}(x, row.names=NULL, optional=TRUE, 
        type=c("response", "link"), ...)
\method{as.data.frame}{efflist}(x, row.names=NULL, optional=TRUE, type, ...)
\method{as.data.frame}{effpoly}(x, row.names=NULL, optional=TRUE, ...)
\method{as.data.frame}{efflatent}(x, row.names=NULL, optional=TRUE, ...)
\method{vcov}{eff}(object, ...)
}

\arguments{
  \item{x, object}{an object consisting of fitted values and other information needed to draw effects plots that is produced by functions in the \code{effects} package.}

  \item{type}{fitted values are by default printed by these functions in the \code{"response"} scale.  For models with a link function like a GLM, fitted values in the linear predictor scale are obtained by setting \code{type="link"}.  For polytomous response models setting \code{type="logits"} returns fitted values in the logit scale.}

  \item{row.names, optional}{arguments to \code{as.data.frame} not used by these methods.}
  
  \item{...}{other arguments passed on}
}

\value{
The \code{print} methods return the fitted values in tables.  The \code{summary} methods return the fitted values and 95 percent condifence intervals, also in tables.  The \code{as.data.frame} method returns fitted values, standard errors, and 95 percent confidence intervals as a data frame, or as a list of data frames for the \code{efflist} method.  The \code{vcov} method returns the covariance matrix of the fitted values.
}


\author{John Fox \email{jfox@mcmaster.ca} and Jangman Hong.}

\examples{

mod.cowles <- glm(volunteer ~ sex + neuroticism*extraversion,
                  data=Cowles, family=binomial)
eff.cowles <- predictorEffects(mod.cowles)
print(eff.cowles)
print(eff.cowles[["neuroticism"]], type="link")
summary(eff.cowles[["neuroticism"]], type="link")
as.data.frame(eff.cowles)
# covariance matrix of fitted values in linear predictor scale
vcov(eff.cowles[[1]]) 
}

\keyword{hplot}
\keyword{models}
