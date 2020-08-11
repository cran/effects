## ----setopts,echo=FALSE-------------------------------------------------------
library("knitr")
opts_chunk$set(fig.width=5,fig.height=5,tidy=TRUE,
               out.width="0.8\\textwidth",echo=TRUE)
options(prompt=" ")

## ----echo=FALSE, results='hide', include=FALSE--------------------------------
#options(continue="+    ", prompt="R> ", width=76)
options(show.signif.stars=FALSE)
options(scipen=3)
library(effects)

## ----include=FALSE------------------------------------------------------------
library(knitr)
opts_chunk$set(
tidy=FALSE,fig.width=5,fig.height=5,cache=FALSE,comment=NA, prompt=TRUE
)
render_sweave()

## ----echo=FALSE, results='hide', include=FALSE----------------------------
options(continue="    ", prompt=" ", width=76)
options(show.signif.stars=FALSE)
options(scipen=3)

## ----fig.height=4,fig.width=8---------------------------------------------
library(effects)
g1 <- lm(prestige ~ education + type + education:type, data = Prestige)
plot(predictorEffects(g1), lines=list(multiline=TRUE))

## -------------------------------------------------------------------------
data(Orthodont, package="nlme")
g2 <- lme4::lmer(distance ~ age + Sex + (1 |Subject), data = Orthodont)
g2

## ----fig.height=4,fig.width=8---------------------------------------------
plot(predictorEffects(g2))

