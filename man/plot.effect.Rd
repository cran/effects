\name{plot.effects}
\alias{plot.effect}
\alias{plot.effects}
\alias{plot.predictoreff}
\alias{plot.predictorefflist}
\alias{plot.eff}
\alias{plot.effpoly}
\alias{plot.efflist}
\alias{plot.mlm.efflist}
\alias{[.efflist}

\title{Plots of Effects and Predictor Effects}

\description{
\code{plot} methods for \code{predictoreff}, \code{predictorefflist}, \code{eff}, \code{efflist} and \code{effpoly} objects created by calls other methods in the \code{effects} package.  The plot arguments were substantially changed in mid-2017. For more details and many examples, see the \href{../doc/predictor-effects-gallery.pdf}{Predictor Effects Graphics Gallery} vignette.
}

\usage{

\method{plot}{eff}(x, x.var, 
    main=paste(effect, "effect plot"),
    symbols=TRUE, lines=TRUE, axes, confint,
    partial.residuals, id, lattice, ...,
    # legacy arguments:
    multiline, z.var, rug, xlab, ylab, colors, cex, lty, lwd,
    ylim, xlim, factor.names, ci.style,
    band.transparency, band.colors, type, ticks,
    alternating, rotx, roty, grid, layout,
    rescale.axis, transform.x, ticks.x, show.strip.values,
    key.args, use.splines,
    residuals.color, residuals.pch, residuals.cex, smooth.residuals,
    residuals.smooth.color, show.fitted, span)
    
\method{plot}{efflist}(x, selection, rows, cols, ask=FALSE, graphics=TRUE, lattice, ...)

\method{plot}{predictoreff}(x, x.var, 
    main = paste(names(x$variables)[1], "predictor effect plot"), ...)

\method{plot}{predictorefflist}(x, selection, rows, cols, ask = FALSE, 
    graphics = TRUE, lattice, ...)

\method{plot}{effpoly}(x, x.var=which.max(levels),
    main=paste(effect, "effect plot"),
    symbols=TRUE, lines=TRUE, axes, confint, lattice, ...,
    # legacy arguments:
    type, multiline, rug, xlab, ylab, colors, cex, lty, lwd,
    factor.names, show.strip.values,
    ci.style, band.colors, band.transparency, style,
    transform.x, ticks.x, xlim,
    ticks, ylim, rotx, roty, alternating, grid,
    layout, key.args, use.splines)
    
\method{plot}{mlm.efflist}(x, ...)
}

\arguments{
  \item{x}{an object of class \code{"predictoreff"}, \code{"predictorefflist"}, \code{"eff"}, \code{"effpoly"}, \code{"efflist"}, \code{"mlm.efflist"}, 	or \code{"summary.eff"}, as appropriate.}

  \item{x.var}{the index (number) or quoted name of the covariate or factor to place on the horizontal axis of each panel of the effect plot. The default is the predictor with the largest number of levels or values.  This argument is ignored with \code{predictoreff} objects.}
  
  \item{main}{the title for the plot, printed at the top; the default title is constructed from the name of the effect.}
  
  \item{symbols}{\code{TRUE}, \code{FALSE}, or an optional list of specifications for plotting symbols; if not given, symbol properties are taken from \code{superpose.symbol} in the lattice theme. See Detailed Argument Descriptions under Details for more information.}
  
  \item{lines}{\code{TRUE}, \code{FALSE}, or an optional list of specifications for plotting lines (and possibly areas); if not given, line properties are taken from \code{superpose.line} in the lattice theme. See Detailed Argument Descriptions under Details for more information.}
  
  \item{axes}{an optional list of specifications for the x and y axes; if not given, axis properties take generally reasonable default values. See Details for more information.}
  
  \item{confint}{an optional list of specifications for plotting confidence regions and intervals; if not given, generally reasonable default values are used. See Detailed Argument Descriptions under Details for more information.}
  
  \item{partial.residuals}{an optional list of specifications for plotting partial residuals for linear and generalized linear models; if not given, generally reasonable default values are used. See Detailed Argument Descriptions under Details for more information, along with the \href{../doc/partial-residuals.pdf}{Effect Displays with Partial Residuals} vignette.}
  \item{id}{an optional list of specifications for identifying points when partial residuals are plotted; if not specified, no points are labelled. See Detailed Argument Descriptions under Details for more information.}
  
  \item{lattice}{an optional list of specifications for various lattice properties, such as legend placement; if not given, generally reasonable default values are used. See Detailed Argument Descriptions under Details for more information.}
  
  \item{selection}{the optional index (number) or quoted name of the effect in an efflist object to be plotted; if not supplied, a menu of high-order terms is presented or all effects are plotted.}
  
  \item{rows, cols}{Number of rows and columns in the ``meta-array'' of plots produced for an \code{efflist} object; if either argument is missing, then the meta-layout will be computed by the \code{plot} method.}
  
  \item{ask}{if \code{selection} is not supplied and \code{ask} is \code{TRUE}, a menu of high-order terms is presented; if \code{ask} is \code{FALSE} (the default), effects for all high-order terms are plotted in an array.}
  
  \item{graphics}{if \code{TRUE} (the default), then the menu of terms to plot is presented in a dialog box rather than as a text menu.}
  
  \item{...}{arguments to be passed down. For \code{"predictoreff"} or \code{"predictorefflist"} objects, the arguments passed down can include all the arguments for \code{"eff"}.}
  
\item{multiline, z.var, rug, xlab, ylab, colors, cex, lty, lwd,
     ylim, xlim, factor.names, ci.style,
     band.transparency, band.colors, ticks,
     alternating, rotx, roty, grid, layout,
     rescale.axis, transform.x, ticks.x, show.strip.values,
     key.args, use.splines, type,
     residuals.color, residuals.pch, residuals.cex, smooth.residuals,
     residuals.smooth.color, show.fitted, span, style}{legacy arguments retained for backwards compatibility; if specified, these will take precedence over the newer list-style arguments described above. See \code{\link{LegacyArguments}} for details.}
}

\details{
Effects plots and predictor effects plots are produced by these methods.  The plots are highly customizable using the optional arguments described here.  For example, effects in a GLM are plotted on the scale of the linear predictor, but the vertical axis is labelled on the response scale. This preserves the linear structure of the model while permitting interpretation on what is usually a more familiar scale. This approach may also be used with linear models, for example to display effects on the scale of the response even if the data are analyzed on a transformed scale, such as log or square-root.  See the \code{axes} argument details below to change the scale to response scale, or to linear predictor scale with tick marks labeled in response scale.

When a factor is on the x-axis, the \code{plot} method for \code{eff} objects connects the points representing the effect by line segments, creating a response ``profile.'' If you wish to suppress these lines, add \code{lty=0} to the \code{lines} argument to the call to \code{plot} (see below and the examples).

In a polytomous multinomial or proportional-odds logit model, by default effects are plotted on the probability scale; they may alternatively be plotted on the scale of the individual-level logits.

All of the arguments to plot objects created by \code{Effect} or \code{allEffects} can also be used with objects created by \code{predictorEffect} or \code{predictorEffects}.  

\bold{Detailed Argument Descriptions}

For more information about these arguments and many examples, see the \href{../doc/predictor-effects-gallery.pdf}{Predictor Effects Graphics Gallery} vignette.

Maximizing the flexibility of these plot commands requires inclusion of a myriad of options.  In an attempt to simplify the use of these options, they have been organized into just a few arguments that each accept a list of specifications as an argument.  In a few cases the named entries in the list are themselves lists.

  Each of the following arguments takes an optional list of specifications; any specification absent from the list
assumes its default value. Some of the list elements are themselves lists, so in complex cases, the argument can take
the form of nested lists.  All of these arguments can also be used on objects created with \code{\link{predictorEffects}}.

\describe{
\item{\code{symbols}}{\code{TRUE}, \code{FALSE}, or a list of options that controls the plotting symbols and their sizes for use with factors;
if \code{FALSE} symbols are suppressed; if \code{TRUE} default values are used:
    \describe{
    \item{\code{pch}}{ploting symbols, a vector of plotting characters, with the default taken from \code{trellis.par.get("superpose.symbol")$pch}, typically a vector of 1s (circles).}
    \item{\code{cex}}{plotting character sizes, a vector of values, with the default taken from \code{trellis.par.get("superpose.symbol")$cex}, typically a vector of 0.8s.}
    }
}
\item{\code{lines}}{\code{TRUE}, \code{FALSE}, or a list that controls the characteristics of lines drawn on a plot, and also whether or not multiple lines should be drawn in the same panel in the plot; if \code{FALSE} lines are suppressed; if \code{TRUE} default values are used:
    \describe{
    \item{\code{multiline}}{display a multiline plot in each panel; the default is \code{TRUE} if there are no standard errors
    in the \code{"eff"} object, \code{FALSE} otherwise. For an \code{"effpoly"} object \code{multline=TRUE} causes all of the response
    levels to be shown in the same panel rather than in separate panels.}
  \item{z.var}{for linear, generalized linear or mixed models,
  	the index (number) or quoted name of the covariate or factor for which
    individual lines are to be drawn in each panel of the effect plot. 
    The default is the
    predictor with the smallest number of levels or values. This argument is only
    used for multipline plots.}
    \item{\code{lty}}{vector of line types, with the default taken from \code{trellis.par.get("superpose.line")$lty}, typically a vector of 1s (solid lines).}
    \item{\code{lwd}}{vector of line widths, with the default taken from \code{trellis.par.get("superpose.line")$lwd}, typically a vector with 2 in the first position followed by 1s.}
    \item{\code{col}}{a vector of line colors, with the default taken from from \code{trellis.par.get("superpose.line")$col}, used both for lines and
    for areas in stacked area plots for \code{"effpoly"} objects; in the latter case, the default colors for an ordered response are instead generated by
    \code{\link[colorspace]{sequential_hcl}} in the \pkg{colorspace} package.}
    \item{\code{splines}}{use splines to smooth plotted effect lines; the default is \code{TRUE}.}
    }
}
\item{\code{axes}}{a list with elements \code{x}, \code{y}, \code{alternating}, and \code{grid} that control axis limits, ticks, and labels.
The \code{x} and \code{y} elements may themselves be lists.

The \code{x} entry is a list with elements named for predictors, with each predictor element itself a list with the following elements:
    \describe{
    \item{\code{lab}}{axis label, defaults to the name of the predictor; may either
    be a text string or a list with the text label (optionally named \code{label}) 
    as its first element and the named element \code{cex} as its second element.}
    \item{\code{lim}}{a two-element vector giving the axis limits, with the default determined from the data.}
    \item{\code{ticks}}{a list with either element \code{at}, a vector specifying locations for the ticks marks, or \code{n}, the number
        of tick marks.}
    \item{\code{transform}}{transformations to be applied to the horizontal axis of a numeric predictor,
        in the form of a list of two functions, with
        element names \code{trans} and \code{inverse}. The
        \code{trans} function is applied to the values of the predictor, and \code{inverse} is used for computing
        proper axis tick labels. The default is not to transform the predictor axis.}
    }

Two additional elements may appear in the \code{x} list, and apply to all predictors:
    \describe{
    \item{\code{rotate}}{angle in degrees to rotate tick labels; the default is 0.}
    \item{\code{rug}}{display a rug plot showing the marginal distribution of a numeric predictor; the default is \code{TRUE}.}
    }

The \code{y} list contains \code{lab}, \code{lim}, \code{ticks}, and \code{rotate} elements
(similar to those specified for individual predictors in the \code{x} list), along with the additional \code{type}, \code{transform}, and \code{style} elements:
    \describe{
    \item{\code{type}}{for plotting linear or generalized linear models, \code{"rescale"} (the default) plots the vertical
    axis on the link scale (e.g., the logit scale for a logit model) but labels the axis on the response
    scale (e.g., the probability scale for a logit model);
    \code{"response"} plots and labels the vertical axis on the scale of the response (e.g., the probability scale for a logit model); and
    \code{"link"} plots and labels the vertical axis on the scale of the link (e.g., the logit scale for a logit model).
    For polytomous logit models, this element is either \code{"probability"} or \code{"logit"}, with the former as the default.}
    \item{\code{transform}}{primarily for linear or linear mixed models, this argument is used to apply an arbitrary transformation to the vertical axis.  For example, if fitting a linear model with response \code{log(y)}, then setting \code{transform=exp} would plot \code{exp(log(y)) = y} on the vertical axis.  If the response were \code{1/y}, then use \code{transform=function(yt) 1/yt}, since the reciprocal is its own inverse.  The \code{transform} argument can also be a list of two functions.  For example with a response \code{log(y)}, the specification \code{transform=list(trans=log, inverse=log), type="rescale"} will plot in log-scale, but will label tick marks in arithmetic scale; see the example below.  The specification \code{transform=list(trans=log, inverse=exp), type="response"} is equivalent to \code{transform=exp}.   When \code{type="response"} the \code{lab} argument will geneally be used to get a label for the axis that matches the untransformed response.  If this argument is used with a generalized linear model or another model with a non-identity link function, the function is applied to the linear predictor, and will probably not be of interest.}
    \item{\code{style}}{for polytomous logit models, this element can take on the value \code{"lines"}
    (the default) or \code{"stacked"} for line plots or stacked-area plots, respectively.}
    }

Other elements:
    \describe{
    \item{\code{alternating}}{if \code{TRUE} (the default), the tick labels alternate by panels in
    multi-panel displays from left to right and top to bottom; if \code{FALSE}, tick labels
    appear at the bottom and on the left.}
    \item{\code{grid}}{if \code{TRUE} (the default is \code{FALSE}), add grid lines to the plot.}
    }
}

\item{\code{confint}}{specifications to add/remove confidence intervals or regions from a plot, and to set the nominal confidence level.

    \describe{
    \item{\code{style}}{one of \code{"auto"}, \code{"bars"}, \code{"lines"}, \code{"bands"}, and \code{"none"}; the default
    is \code{"bars"} for factors, \code{"bands"} for numeric predictors, and \code{"none"} for multiline plots; \code{"auto"} also produces \code{"bars"} for factors
    and \code{"bands"} for numeric predictors, even in multiline plots.}
    \item{\code{alpha}}{transparency of confidence bands; the default is 0.15.}
    \item{\code{col}}{colors; the default is taken from the line colors.}
    }
}

\item{\code{partial.residuals}}{specifications concerning the addition of partial residuals to the plot.

    \describe{
    \item{\code{plot}}{display the partial residuals;
    the default is \code{TRUE} if residuals are present in the \code{"eff"} object, \code{FALSE} otherwise.}
    \item{\code{fitted}}{show fitted values as well as residuals; the default is \code{FALSE}.}
    \item{\code{col}}{color for partial residuals; the default is the second line color.}
    \item{\code{pch}}{plotting symbols for partial residuals; the default is 1, a circle.}
    \item{\code{cex}}{size of symbols for partial residuals; the default is 1.}
    \item{\code{smooth}}{draw a loess smooth of the partial residuals; the default is \code{TRUE}.}
    \item{\code{span}}{span for the loess smooth; the default is 2/3.}
    \item{\code{smooth.col}}{color for the loess smooth; the default is the second line color.}
    \item{\code{lty}}{line type for the loess smooth; the default is the first line type, normally 1 (a solid line).}
    \item{\code{lwd}}{line width for the loess smooth; the default is the first line width, normally 2.}
    }

}

\item{\code{id}}{specifications for optional point identification when partial residuals are plotted.

  \describe{
  \item{\code{n}}{number of points to identify; default is \code{2} if \code{id=TRUE}
    and \code{0} if \code{id=FALSE}. Points are selected based on the Mahalanobis
    distances of the pairs of x-values and partial residuals from their centroid.}
  \item{\code{col}}{color for the point labels; default is the same as the color of the partial residuals.}
  \item{\code{cex}}{relative size of text for point labels; default is \code{0.75}.}
  \item{\code{labels}}{vector of point labels; the default is the names of the residual vector,
  which is typically the row names of the data frame to which the model is fit.}
  }

}

\item{\code{lattice}}{the plots are drawn with the \pkg{\link{lattice}} package, generally by the \code{\link{xyplot}} function.  These specifications are passed as arguments to the functions that actually draw the plots.

    \describe{
    \item{\code{layout}}{the \code{layout} argument to the \pkg{lattice} function \code{\link{xyplot}}
  	(or, in some cases \code{\link{densityplot}}), which
  	is used to draw the effect display; if not specified, the plot will be formatted so that
  	it appears on a single page.}
  	\item{\code{key.args}}{a key, or legend, is added to the plot if \code{multiline=TRUE}.  This argument is a list with components that determine the the placement and other characteristics of the key.  The default if not set by the user is \code{key.args = list(space="top", columns=2, border=FALSE, fontfamily="serif", cex.title=.80, cex=0.75)}.  If there are more than 6 groups in the plot, \code{columns} is set to 3. For stacked-area plots, the default is a one-column key.  In addition to the arguments shown explicitly below, any of the arguments listed in the \code{\link[lattice]{xyplot}} documentation in the \code{key} section can be used.
  	\describe{
  	\item{\code{space}}{determines the placement of the key outside the plotting area, with default \code{space="above"} for above the plot and below its title.  Setting \code{space="right"} uses space to the right of the plot for the key.}
  	\item{\code{x, y, corner}}{used to put the key on the graph itself.  For example, \code{x=.05, y=.95, corner=c(0,1)} will locate the upper-left corner of the key at (.05, .95), thinking of the graph as a unit square.}
  	\item{\code{columns}}{number of columns in the key.  If \code{space="top"}, columns should be 2, 3 or 4; if \code{space="right"}, set \code{columns=1}.}
  	\item{\code{border}}{if \code{TRUE} draw a border around the key; omit the border if \code{FALSE}.}
  	\item{\code{fontfamily}}{the default is \code{"sans"} for the sans-serif font used in the rest of the plot; the alternative is \code{"serif"} for a serif font.}
  	\item{\code{cex, cex.title}}{the default relative size of the font for labels and the title, respectively.  To save space set these to be smaller than 1.}
  	}
}
    \item{\code{strip}}{a list with three elements: \code{factor.names}, which if \code{TRUE}, the default, shows conditioning
    variable names in the panel headers;  \code{values}, which if \code{TRUE}, the default unless partial residuals are plotted,
    displays conditioning variable values in the panel headers, and \code{cex}, the relative size of the text displayed in the strip.}
    \item{\code{array}}{a list with elements \code{row}, \code{col}, \code{nrow}, \code{ncol}, and \code{more},
    used to graph an effect as part of an array of plots; \code{row}, \code{col}, \code{nrow}, and \code{ncol} are used to compose
    the \code{split} argument and \code{more} the \code{more} argument to \code{\link{print.trellis}}.
    The \code{array} argument is automatically set by \code{plot.efflist} and will be ignored if used with that function.}
}

}
}
}

\value{
  The \code{summary} method for \code{"eff"} objects returns a \code{"summary.eff"} object with the following components
  (those pertaining to confidence limits need not be present):
  \item{header}{a character string to label the effect.}
  \item{effect}{an array containing the estimated effect.}
  \item{lower.header}{a character string to label the lower confidence limits.}
  \item{lower}{an array containing the lower confidence limits.}
  \item{upper.header}{a character string to label the upper confidence limits.}
  \item{upper}{an array containing the upper confidence limits.}

  The \code{plot} method for \code{"eff"} objects returns a \code{"plot.eff"} object (an enhanced \code{"trellis"} object); the provided
\code{\link{print}} method plots the object.

  The \code{[} method for \code{"efflist"} objects is used to subset an \code{"efflist"} object and returns an object of the same class.
}


\author{John Fox \email{jfox@mcmaster.ca} and Jangman Hong.}


\seealso{\code{\link{LegacyArguments}}, \code{\link{effect}}, \code{\link{allEffects}}, \code{\link{effectsTheme}},
    \code{\link{xyplot}}, \code{\link{densityplot}}, \code{\link{print.trellis}}, \code{\link{loess}},
	\code{\link[colorspace]{sequential_hcl}},
	and the \href{../doc/predictor-effects-gallery.pdf}{Predictor Effects Graphics Gallery} and \href{../doc/partial-residuals.pdf}{Effect Displays with Partial Residuals} vignettes.
}

\examples{
# also see examples in ?effect

# plot predictorEffects
mod <- lm(prestige ~ education + log(income)*type + women, Prestige)
plot(predictorEffects(mod, ~ income), axes=list(grid=TRUE))
plot(predictorEffects(mod, ~ income), lines=list(multiline=TRUE),
                                                 axes=list(grid=TRUE))
plot(predictorEffects(mod, ~ type), lines=list(multiline=TRUE),
                                                 axes=list(grid=TRUE),
                                                 confint=list(style="bars"))

mod.cowles <- glm(volunteer ~ sex + neuroticism*extraversion,
                  data=Cowles, family=binomial)
eff.cowles <- allEffects(mod.cowles, xlevels=list(extraversion=seq(0, 24, 6)))
eff.cowles
as.data.frame(eff.cowles[[2]]) # neuroticism*extraversion interaction

plot(eff.cowles, 'sex', axes=list(grid=TRUE,
                                  y=list(lab="Prob(Volunteer)"),
                                  x=list(rotate=90)),
                        lines=list(lty=0))

plot(eff.cowles, 'neuroticism:extraversion',
     axes=list(y=list(lab="Prob(Volunteer)",
        ticks=list(at=c(.1,.25,.5,.75,.9)))))

plot(Effect(c("neuroticism", "extraversion"), mod.cowles,
            se=list(type="Scheffe"),
            xlevels=list(extraversion=seq(0, 24, 6))),
     axes=list(y=list(lab="Prob(Volunteer)",
        ticks=list(at=c(.1,.25,.5,.75,.9)))))

\donttest{
    # change color of the confidence bands to 'black' with .15 transparency
plot(eff.cowles, 'neuroticism:extraversion',
     axes=list(y=list(lab="Prob(Volunteer)",
                      ticks=list(at=c(.1,.25,.5,.75,.9)))),
     confint=list(col="red", alpha=.3))

plot(eff.cowles, 'neuroticism:extraversion',
     lines=list(multiline=TRUE),
     axes=list(y=list(lab="Prob(Volunteer)")),
     lattice=list(key.args = list(x = 0.65, y = 0.99, corner = c(0, 1))))

# use probability scale in place of logit scale, all lines are black.
plot(eff.cowles, 'neuroticism:extraversion',
     lines=list(multiline=TRUE, lty=1:8, col="black"),
     axes=list(y=list(type="response", lab="Prob(Volunteer)")),
     lattice=list(key.args = list(x = 0.65, y = 0.99, corner = c(0, 1))),
     confint=list(style="bands"))

plot(effect('sex:neuroticism:extraversion', mod.cowles,
            xlevels=list(extraversion=seq(0, 24, 6))),
     lines=list(multiline=TRUE))

plot(effect('sex:neuroticism:extraversion', mod.cowles,
            xlevels=list(extraversion=seq(0, 24, 6))),
     lines=list(multiline=TRUE),
     axes=list(y=list(type="response")),
     confint=list(style="bands"),
     lattice=list(key.args = list(x=0.75, y=0.75, corner=c(0, 0))))
}

if (require(nnet)){
    mod.beps <- multinom(vote ~ age + gender + economic.cond.national +
                             economic.cond.household + Blair + Hague + Kennedy +
                             Europe*political.knowledge, data=BEPS)
    \donttest{
    plot(effect("Europe*political.knowledge", mod.beps,
                xlevels=list(political.knowledge=0:3)))
    }

    plot(effect("Europe*political.knowledge", mod.beps,
                xlevels=list(political.knowledge=0:3),
                fixed.predictors=list(given.values=c(gendermale=0.5))),
         axes=list(y=list(style="stacked"), x=list(rug=FALSE), grid=TRUE),
         lines=list(col=c("blue", "red", "orange")))
}

if (require(MASS)){
    mod.wvs <- polr(poverty ~ gender + religion + degree + country*poly(age,3),
                    data=WVS)
    plot(effect("country*poly(age, 3)", mod.wvs))

    \donttest{
    plot(effect("country*poly(age, 3)", mod.wvs), lines=list(multiline=TRUE))
    plot(effect("country*poly(age, 3)", mod.wvs),
         axes=list(y=list(style="stacked")),
         lines=list(col=c("gray75", "gray50", "gray25")))

    plot(effect("country*poly(age, 3)", latent=TRUE, mod.wvs))
    }
}

mod.pres <- lm(prestige ~ log(income, 10) + poly(education, 3) + poly(women, 2),
               data=Prestige)
eff.pres <- allEffects(mod.pres)
\donttest{
plot(eff.pres)
plot(eff.pres[1:2])
}
plot(eff.pres[1],
     axes=list(x=list(income=list(transform=list(
         trans=log10, inverse=function(x) 10^x),
         ticks=list(at=c(1000, 2000, 5000, 10000, 20000))))))
         
mod <- lm(log(prestige) ~ income:type + education, data=Prestige)
p1 <- predictorEffects(mod, ~ income)
# log-scale for response
plot(p1, lines=list(multiline=TRUE)) 
# log-scale, with arithmetic tick marks
plot(p1, lines=list(multiline=TRUE), 
     axes=list(y=list(transform=list(trans=log, inverse = exp), 
                      lab="prestige", type="rescale")))
# arithmetic scale and tick marks, with other arguments
plot(p1, lines=list(multiline=TRUE), grid=TRUE,
     lattice=list(key.args=list(space="right", border=TRUE)),
     axes=list(y=list(transform=exp, lab="prestige")))
}

\keyword{hplot}
\keyword{models}
