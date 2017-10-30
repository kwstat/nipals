\name{uscrime}
\alias{uscrime}
\docType{data}
\title{
U.S. Crime rates per 100,00 people
}
\description{
U.S. Crime rates per 100,00 people for 7 categories in each of the 50 U.S. states in 1977.
}
\usage{data("uscrime")}
\format{
  A data frame with 50 observations on the following 8 variables.
  \describe{
    \item{\code{state}}{U.S. state}
    \item{\code{murder}}{murders}
    \item{\code{rape}}{rapte}
    \item{\code{robbery}}{robbery}
    \item{\code{assault}}{assault}
    \item{\code{burglary}}{burglary}
    \item{\code{larceny}}{larceny}
    \item{\code{autotheft}}{automobile thefts}
  }
}
\details{
Two missing values.
}
\source{
http://documentation.sas.com/api/docsets/stathpug/14.2/content/stathpug_code_hppriex3.htm?locale=en
}
\references{
  None
}
\examples{
  dat <- uscrime[complete.cases(uscrime), ]
  dat <- as.matrix(dat[ , -1])
  m1 <- nipals(dat) # complete-data method
}
\keyword{datasets}
