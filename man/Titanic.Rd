\name{TitanicSurvival}
\alias{TitanicSurvival}
\docType{data}
\title{Survival of Passengers on the Titanic}
\description{
  Information on the survival status, sex, age, and passenger class of 1309
  passengers in the Titanic disaster of 1912.  
}

\usage{TitanicSurvival}

\format{
  A data frame with 1309 observations on the following 4 variables.
  \describe{
    \item{\code{survived}}{\code{no} or \code{yes}.}
    \item{\code{sex}}{\code{female} or \code{male}}
    \item{\code{age}}{in years (and for some children, fractions of a year); age
      is missing for 263 of the passengers.}
    \item{\code{passengerClass}}{\code{1st}, \code{2nd}, or \code{3rd} class.}
  }
}

\details{
  This is part of a larger data set compiled by Thomas Cason. 
  Many additional details are given in the sources cited below.
}

\source{
  Data set \code{titanic3} from 
  \url{http://biostat.mc.vanderbilt.edu/twiki/bin/view/Main/DataSets}.
}

\references{
  \url{http://www.encyclopedia-titanica.org/}
  
  F. E. Harrell, Jr. (2001)
  \emph{Regression Modeling Strategies}
  New York: Springer.
}
\examples{
summary(TitanicSurvival)
  
titanic <- glm(survived ~ (passengerClass + sex + age)^2, 
                data=TitanicSurvival, family=binomial)
  
titanic.all <- allEffects(titanic, typical=median, 
	given.values=c(passengerClass2nd=1/3, passengerClass3rd=1/3, sexmale=0.5))
  	
plot(titanic.all, ticks=list(at=c(.01, .05, seq(.1, .9, by=.2), .95, .99)), 
    ask=FALSE)

plot(effect("passengerClass*sex*age", titanic, xlevels=list(age=0:65)), 
	ticks=list(at=c(.001, .005, .01, .05, seq(.1, .9, by=.2), .95, .99, .995))) 
}
\keyword{datasets}
