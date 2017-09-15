\name{effectsTheme}
\alias{effectsTheme}

\title{Set the lattice Theme for Effect Plots}

\description{
Set the \pkg{lattice} theme (see \code{\link[lattice]{trellis.device}}) appropriately for effect plots. This function is invoked
automatically when the \pkg{effects} package is loaded \emph{if} the \pkg{lattice} package hasn't previously been loaded. A typical
call is \code{lattice::trellis.par.set(effectsTheme())}.
}

\usage{
effectsTheme(strip.background = list(col = gray(seq(0.95, 0.5, length = 3))), 
    strip.shingle = list(col = "black"), clip = list(strip = "off"), 
    superpose.line = list(lwd = c(2, rep(1, 6))))
}

\arguments{
  \item{strip.background}{colors for the background of conditioning strips at the top of each panel; the default uses shades of gray and 
    makes allowance for up to three conditioning variables.}
  \item{strip.shingle}{when lines rather than numeric values are used to indicate the values of conditioning variables, the default sets the color of
    the lines to black.}
  \item{clip}{the default allows lines showing values of conditioning variables to extend slightly beyond the boundaries of 
    the strips---making the lines more visible at the extremes.}
  \item{superpose.line}{the default sets the line width of the first (of seven) lines to 2.}
}

\value{
a list suitable as an argument for \code{\link[lattice]{trellis.par.set}}; current values of modified parameters are
supplied as an attribute.
}

\author{John Fox \email{jfox@mcmaster.ca}}

\seealso{\code{\link[lattice]{trellis.device}}, \code{\link[lattice]{trellis.par.set}}}

\examples{
\dontrun{
lattice::trellis.par.set(effectsTheme())
}
}


\keyword{utilities}
\keyword{device}
