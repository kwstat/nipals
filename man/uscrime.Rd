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
  There are two missing values.
}
\source{
  Documentation Example 3 for PROC HPPRINCOMP.
  http://documentation.sas.com/api/docsets/stathpug/14.2/content/stathpug_code_hppriex3.htm?locale=en
}
\references{
  SAS/STAT User's Guide: High-Performance Procedures.
  The HPPRINCOMP Procedure.
  http://support.sas.com/documentation/cdl/en/stathpug/67524/HTML/default/viewer.htm#stathpug_hpprincomp_toc.htm
}

\examples{
  library(nipals)
  data(uscrime)

  # SAS deletes rows with missing values
  dat <- uscrime[complete.cases(uscrime), ]
  dat <- as.matrix(dat[ , -1])
  m1 <- nipals(dat) # complete-data method

  # Traditional NIPALS with missing data  
  dat <- uscrime
  dat <- as.matrix(dat[ , -1])
  m2 <- nipals(dat, gramschmidt=FALSE) # missing 
  round(crossprod(m2$loadings),3) # Prin Comps not quite orthogonal
  
  # Gram-Schmidt corrected NIPALS
  m3 <- nipals(dat, gramschmidt=TRUE) # TRUE is default
  round(crossprod(m3$loadings),3) # Prin Comps are orthogonal
}
\keyword{datasets}
