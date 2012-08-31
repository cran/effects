\name{WVS}
\alias{WVS}
\docType{data}
\title{World Values Surveys}
\description{
  Data from the World Values Surveys 1995-1997 for Australia, Norway, Sweden, and the United States.
}
\usage{WVS}
\format{
  A data frame with 5381 observations on the following 6 variables.
  \describe{
    \item{\code{poverty}}{``Do you think that what the government is doing for people in poverty in
	this country is about the right amount, too much, or too little?'' (ordered): 
	\code{Too Little}, \code{About Right}, \code{Too Much}}.
    \item{\code{religion}}{Member of a religion: \code{no} or \code{yes}.}
    \item{\code{degree}}{Held a university degree: \code{no} or \code{yes}.}
    \item{\code{country}}{\code{Australia}, \code{Norway}, \code{Sweden}, or \code{USA}.}
    \item{\code{age}}{in years.}
    \item{\code{gender}}{\code{male} or \code{female}.}
  }
}

\references{
  J. Fox and R. Andersen (2006)
   Effect displays  for multinomial and proportional-odds logit models.
  \emph{Sociological Methodology} \bold{36}, 225--255.
}

\examples{
summary(WVS)

\donttest{ 
require(splines)  # for bs()
wvs <- polr(poverty ~ gender + country*(religion + degree + bs(age, 4)), data=WVS)

plot(effect("country*bs(age,4)", wvs, xlevels=list(age=18:83),
  given.values=c(gendermale=0.5)), rug=FALSE) 
plot(effect("country*bs(age,4)", wvs, xlevels=list(age=18:83),
  given.values=c(gendermale=0.5)), rug=FALSE, style="stacked")
plot(effect("country*bs(age,4)", wvs, xlevels=list(age=18:83),
  given.values=c(gendermale=0.5), latent=TRUE), rug=FALSE)
}
}

\keyword{datasets}
