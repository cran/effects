\name{EffectMethods}
\alias{Effect.default}
\alias{Effect.merMod}
\alias{Effect.rlmerMod}
\alias{Effect.lme}
\alias{Effect.gls}
\alias{Effect.poLCA}
\alias{Effect.clm2}
\alias{Effect.clm}
\alias{Effect.clmm}
\alias{Effect.mlm}
\alias{Effect.betareg}

\title{Functions For Constructing Effect Displays for Many Modeling Paradigms}
\description{The \code{Effect}, \code{effect} and \code{predictorEffects} methods are used to draw effects plots to visualize a fitted regression surface.  These plots can be drawn at least in principle for any model that uses a linear predictor.  Methods for modeling paradigms than the basic \code{lm}, \code{glm}, \code{multinom} and \code{polr} methods are documented here. To add methods for new classes of statistical models, see the vignette \href{../doc/adding-Effect-methods.pdf}{Defining Effect Methods for Other Models}.
}

\usage{

\method{Effect}{default}(focal.predictors, mod,  ...,
    sources=NULL)

\method{Effect}{gls}(focal.predictors, mod, ...)

\method{Effect}{clm2}(focal.predictors, mod, ...)

\method{Effect}{clmm}(focal.predictors, mod, ...)

\method{Effect}{clm}(focal.predictors, mod, ...)

\method{Effect}{merMod}(focal.predictors, mod, ...,
    KR=FALSE)

\method{Effect}{rlmerMod}(focal.predictors, mod, ...)

\method{Effect}{lme}(focal.predictors, mod, ...)

\method{Effect}{poLCA}(focal.predictors, mod, ...)

\method{Effect}{mlm}(focal.predictors, mod, response, ...)

\method{Effect}{betareg}(focal.predictors, mod,  ...)
}

\arguments{
  \item{focal.predictors}{a character vector of one or more predictors in the model in any order.}
  \item{mod}{a fitted model object of the appropriate class. }
  \item{...}{additional arguments passed to other \code{Effect}.  See \code{\link{Effect}} for all the arguments included.}
  \item{response}{for an \code{"mlm"} object, a vector containing the name(s) or indices of one or more response variable(s). The default is to use all responses in the model.}
  \item{sources}{This argument appears only in the default method for \code{Effect}, and allows the user to draw effects plots for fitting methods for which there are not existing methods in the effects package.  Seven arguments are provided:

    \describe{
  	  \item{type}{the default is \code{"glm"}, which assumes the modeling method shares characteristics with a generalized linear model, including a univariate response, a linear predictor, and possibly a error family and link function.}

      \item{call}{For S3 objects, the default is \code{ojbect$call}, returning the call that created the object.  This is used to harvest standard arguments like \code{data}, \code{subset} and \code{family}.}

      \item{formula}{the formula for the linear predictor, defaulting to \code{formula(object)}.}

      \item{family}{if the model object includes an error family, but it is not returned by \code{family(object)}, specify the family with this argument; otherwise it can be ignored.}

      \item{method}{For ordinal response models only, see the \code{method} argument to \code{\link{polr}}.}

      \item{coefficients}{The estimates of the coefficents in the linear predictor, with default \code{coef(object)}.}

      \item{vcov}{the estimated variance covariance materix to be used in computing errors in the effects plots; default is code{vcov(object)}.}
    }
  }
\item{KR}{if \code{TRUE} and the \pkg{pbkrtest} package is installed,
    use the Kenward-Roger coefficient covariance
    matrix to compute effect standard errors for linear mixed models fit with \code{\link[lme4]{lmer}}
    in the \pkg{lme4} package.  The default is \code{FALSE} because the computation can be
    very slow.}
  }

\details{
  Most of these methods simply call the \code{Effect.default} method with the appropriate values in the arguement \code{sources}.  See the vignettte Effect Methods in the vignettes for the effects package.  All the iteresting work is done by the methods described in \code{\link{Effect}}.}

\value{See \code{\link{Effect}}}


\references{
    The \href{../doc/adding-Effect-methods.pdf}{Defining Effect Methods for Other Models} vignette.
    }



\author{John Fox \email{jfox@mcmaster.ca}, Sanford Weisberg \email{sandy@umn.edu}}


\seealso{\code{\link{Effect}} and the links therein.}

\examples{
\dontrun{
# lme
require(nlme)
fm1 <- lme(distance ~ age + Sex, data = Orthodont, random = ~ 1)
plot(predictorEffects(fm1))

# gls
library(nlme)
g <- gls(Employed ~ GNP + Population,
         correlation=corAR1(form= ~ Year), data=longley)
print(predictorEffects(g))

# lmer  uses method Effect.lmerMod
if("package:nlme" %in% search()) detach(package:nlme)
require(lme4)
data("Orthodont", package="nlme")
fm2 <- lmer(distance ~ age + Sex + (1 |Subject), data = Orthodont)
plot(allEffects(fm2))

# glmer  uses method Effect.lmerMod
require(lme4)
gm1 <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
              data = cbpp, family = binomial)
as.data.frame(predictorEffect("period", gm1))

# rlmer  uses method Effect.rlmerMod
require(lme4)
fm3 <- robustlmm::rlmer(distance ~ age + Sex + (1 |Subject), data = Orthodont)
plot(effect("age:Sex", fm3))
plot(predictorEffects(fm3, ~ age + Sex))

# betareg from the betareg package
library(betareg)
library(lme4)
data("GasolineYield", package = "betareg")
gy_logit <- betareg(yield ~ batch + temp, data = GasolineYield)
summary(gy_logit)
Effect("batch", gy_logit)
predictorEffects(gy_logit)

# clm in ordinal
require(ordinal)
require(MASS)
mod.wvs1 <- clm(poverty ~ gender + religion + degree + country*poly(age,3),
    data=WVS)
plot(Effect(c("country", "age"), mod.wvs1), 
     lines=list(multiline=TRUE), layout=c(2, 2))

# clm2
require(ordinal)
require(MASS)
v2 <- clm2(poverty ~ gender + religion + degree + country*poly(age,3),data=WVS)
plot(emod2 <- Effect(c("country", "age"), v2))

# clmm
require(ordinal)
require(MASS)
mm1 <- clmm(SURENESS ~ PROD + (1|RESP) + (1|RESP:PROD), 
        data = soup, link = "logit", threshold = "flexible")
plot(Effect("PROD", mm1),lines=list(multiline=TRUE))

# poLCA
library(poLCA)
data(election)
f2a <- cbind(MORALG,CARESG,KNOWG,LEADG,DISHONG,INTELG,
             MORALB,CARESB,KNOWB,LEADB,DISHONB,INTELB)~PARTY
nes2a <- poLCA(f2a,election,nclass=3,nrep=5)    # log-likelihood: -16222.32
allEffects(nes2a)

# multivariate linear model
data(Baumann, package="carData")
b1 <- lm(cbind(post.test.1, post.test.2, post.test.3) ~ group + 
    pretest.1 + pretest.2, data = Baumann))
plot(Effect("group", b1)

}
}

\keyword{hplot}
\keyword{models}
