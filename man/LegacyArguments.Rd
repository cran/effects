\name{LegacyArguments}
\alias{LegacyArguments}
\alias{Legacy Arguments}

\title{Legacy Arguments for \code{plot} and \code{Effect} Methods}

\description{
Prior to verson 4.0-0 of the \pkg{effects} package, there were many (literally dozens) of
arguments to the \code{plot} methods for \code{"eff"} and \code{"effpoly"} objects. 

In version 4.0-0 of the package, we have consolidated these arguments into a much smaller
number of arguments (e.g., \code{lines}, \code{points}, \code{axes}) that take lists of
specifications. We have similarly consolidated some of the arguments to \code{Effect} methods
into the \code{confint} and \code{fixed.predictors} arguments.

For backwards compatibility, we have to the extent possible retained the
older arguments. If specified, these legacy arguments take precedence over the newer
list-style arguments
}

\details{
Here is the correspondence between the old and new arguments.

For \code{plots} methods:

\describe{
\item{\code{multiline=TRUE/FALSE}}{\code{lines=list(multiline=TRUE/FALSE)}}

\item{\code{colors={vector of colors}}}{\code{lines=list(col={vector of colors})}}

\item{\code{lty={vector of line types}}}{\code{lines=list(lty={vector of line types})}}

\item{\code{lwd={vector of line widths}}}{\code{lines=list(lwd={vector of line widths})}}

\item{\code{use.splines=TRUE/FALSE}}{\code{lines=list(splines=TRUE/FALSE)}}

\item{\code{cex={number}}}{\code{points=list(cex={number})}}

\item{\code{rug=TRUE/FALSE}}{\code{axes=list(x=list(rug=TRUE/FALSE)}}

\item{\code{xlab={"axis title"}}}{\code{axes=list(x=list(lab={"axis title"}))}}

\item{\code{xlim={c(min, max)}}}{\code{axes=list(x=list(lim={c(min, max)}))}}

\item{\code{rotx={degrees}}}{\code{axes=list(x=list(rot={degrees}))}}

\item{\code{ticks.x=list({tick specifications})}}{\code{axes=list(x=list(ticks=list({tick specifications})))}}

\item{\code{transform.x=list(link={function}, inverse={function})}}{\code{axes=list(x=list(transform=list({lists of transformations by predictors})))}}

\item{\code{ylab={"axis title"}}}{\code{axes=list(y=list(lab={"axis title"}))}}

\item{\code{ylim={c(min, max)}}}{\code{axes=list(y=list(lim={c(min, max)}))}}

\item{\code{roty={degrees}}}{\code{axes=list(y=list(rot={degrees}))}}

\item{\code{ticks=list({tick specifications})}}{\code{axes=list(y=list(ticks=list({tick specifications})))}}

\item{\code{alternating=TRUE/FALSE}}{\code{axes=list(alternating=TRUE/FALSE)}}

\item{\code{grid=TRUE/FALSE}}{\code{axes=list(grid=TRUE/FALSE)}}

\item{\code{ci.style="bands"/"lines"/"bars"/"none"}}{\code{confint=list(style="bands"/"lines"/"bars"/"none"})}

\item{\code{band.transparency={number}}}{\code{confint=list(alpha={number})}}

\item{\code{band.colors={vector of colors}}}{\code{confint=list(col={vector of colors})}}

\item{\code{residuals.color={color}}}{\code{partial.residuals=list(col={color})}}

\item{\code{residuals.pch={plotting character}}}{\code{partial.residuals=list(pch={plotting character})}}

\item{\code{residuals.cex={number}}}{\code{partial.residuals=list(cex={number})}}

\item{\code{smooth.residuals=TRUE/FALSE}}{\code{partial.residuals=list(smooth=TRUE/FALSE)}}

\item{\code{residuals.smooth.color={color}}}{\code{partial.residuals=list(smooth.col={color})}}

\item{\code{span={number}}}{\code{partial.residuals=list(span={number})}}

\item{\code{show.fitted=TRUE/FALSE}}{\code{partial.residuals=list(fitted=TRUE/FALSE)}}

\item{\code{factor.names=TRUE/FALSE}}{\code{lattice=list(strip=list(factor.names=TRUE/FALSE))}}

\item{\code{show.strip.values=TRUE/FALSE}}{\code{lattice=list(strip=list(values=TRUE/FALSE))}}

\item{\code{layout={lattice layout}}}{\code{lattice=list(layout={lattice layout})}}

\item{\code{key.args={lattice key args}}}{\code{lattice=list(key.args={lattice key args})}}

\item{\code{style="lines"/"stacked"}}{for \code{plot.effpoly}, \code{axes=list(y=list(style="lines"/"stacked"))}}

\item{\code{rescale.axis=TRUE/FALSE}}{\code{type="rescale"/"response"/"link"}}
}

For \code{Effect} methods:

\describe{
\item{\code{confint=TRUE/FALSE} or a list}{may be substituted for the \code{se} argument.}

\item{\code{confidence.level={number}}}{\code{se=list(level={number})}}

\item{\code{given.values={named vector}}}{\code{fixed.predictors=list(given.values={named vector})}}

\item{\code{typical={function}}}{\code{fixed.predictors=list(typical={function})}}

\item{\code{offset={function}}}{\code{fixed.predictors=list(offset={function})}}

\item{\code{partial.residuals=TRUE/FALSE}}{\code{residuals=TRUE/FALSE}}
}

}

\author{John Fox \email{jfox@mcmaster.ca}}


\seealso{
\code{\link{Effect}}, \code{\link{plot.eff}}, \code{\link{plot.effpoly}}
}

\keyword{hplot}
