\name{plot.predictoreff}
\alias{plot.predictoreff}
\alias{plot.predictorefflist}

\title{
Draw Predictor Effect Plots
}
\description{
These functions call \code{\link{plot.eff}} and \code{\link{plot.efflist}} to draw predictor effect plots.
}
\usage{


\method{plot}{predictoreff}(x, x.var, 
    main = paste(names(x$variables)[1], "predictor effect plot"), ...)

\method{plot}{predictorefflist}(x, selection, rows, cols, ask = FALSE, graphics = TRUE, 
     lattice, ...)
}
\arguments{
  \item{x}{
An object of class \code{predictoreff} or \code{predictorefflist}.
}
  \item{x.var}{
the index (number) or quoted name of the covariate or factor to place on the horizontal axis of each panel of the effect plot. The default is the predictor with the largest number of levels or values.
}
  \item{main}{
the title for the plot, printed at the top; the default title is constructed from the name of the effect.
}
  \item{\dots}{
arguments to be passed to \code{\link{plot.eff}} or \code{\link{plot.efflist}}.
}
  \item{selection}{
the optional index (number) or quoted name of the effect in an effect list to be plotted; if not supplied, a menu of high-order terms is presented or all effects are plotted.
}
  \item{rows, cols}{Number of rows and columns in the "meta-array"" of plots produced for an efflist object; if either argument is missing, then the meta-layout will be computed by the plot method.
}
  \item{ask}{if selection is not supplied and ask is \code{TRUE}, a menu of high-order terms is presented; if ask is \code{FALSE} (the default), effects for all high-order terms are plotted in an array.
}
  \item{graphics}{
  if \code{TRUE} (the default), then the menu of terms to plot is presented in a dialog box rather than as a text menu.
  }
  \item{lattice}{
  argument passed to \code{\link{plot.efflist}}.
  }
}
\details{
The \code{plot.predictoreff} calls the method \code{plot.eff} and \code{plot.predictorefflist} calls \code{plot.efflist}.  Both of these functions are documented at \code{\link{plot.eff}}.  Warning:  By default, the functions documented here use the argument \code{lines=list(multiline=TRUE)} while direct calls to the underlying functions use \code{lines=list(multiline = FALSE)} if standard errors were computed by the call to create the object \code{x}.
}
\value{
See the documentation for \code{\link{plot.eff}}.
}

\author{
S. Weisberg, \email{sandy@umn.edu}
}


\seealso{
\code{\link{predictorEffect}}, \code{\link{plot.eff}}.
}
\examples{
mod <- lm(prestige ~ type*(education + income + women), Prestige)
plot(predictorEffects(mod, ~ income))
}

\keyword{hplot}
\keyword{models}
