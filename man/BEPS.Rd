\name{BEPS}
\alias{BEPS}
\docType{data}
\title{British Election Panel Study}
\description{
  These data are drawn from the 1997-2001 British Election Panel Study (BEPS).
}
\usage{BEPS}
\format{
  A data frame with 1525 observations on the following 10 variables.
  \describe{
    \item{\code{vote}}{Party choice: \code{Conservative}, \code{Labour}, or \code{Liberal Democrat}}
    \item{\code{age}}{in years}
    \item{\code{economic.cond.national}}{Assessment of current national economic conditions, 1 to 5.}
    \item{\code{economic.cond.household}}{Assessment of current household economic conditions, 1 to 5.}
    \item{\code{Blair}}{Assessment of the Labour leader, 1 to 5.}
    \item{\code{Hague}}{Assessment of the Conservative leader, 1 to 5.}
    \item{\code{Kennedy}}{Assessment of the leader of the Liberal Democrats, 1 to 5.}
    \item{\code{Europe}}{an 11-point scale that measures respondents' attitudes
		toward European integration. High scores represent `Eurosceptic'
		sentiment.}
    \item{\code{political.knowledge}}{Knowledge of parties' positions on European integration, 0 to 3.}
    \item{\code{gender}}{\code{female} or \code{male}.}
  }
}

\references{
  J. Fox and R. Andersen (2006)
   Effect displays  for multinomial and proportional-odds logit models.
  \emph{Sociological Methodology} \bold{36}, 225--255.
}

\examples{
summary(BEPS)

\donttest{  
require(splines)  # for bs()
beps <- multinom(vote ~ age + gender + economic.cond.national + economic.cond.household 
	+ Blair + Hague + Kennedy + bs(Europe, 3)*political.knowledge, data=BEPS)
europe.knowledge <- effect("bs(Europe, 3)*political.knowledge", beps, 
	xlevels=list(Europe=seq(1, 11, length=50), political.knowledge=0:3),
  	given.values=c(gendermale=0.5))
plot(europe.knowledge)
plot(europe.knowledge, style="stacked", colors=c("blue", "red", "orange"), rug=FALSE)
}
}

\keyword{datasets}
