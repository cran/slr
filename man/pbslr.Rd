\name{pbslr}
\alias{pbslr}
\title{Partially balanced semi-Latin rectangles of cell size two}
\description{This function generates a partially balanced semi-Latin rectangles of cell size two}
\usage{pbslr(v)}
\arguments{
	\item{v}{number of treatments}	
	}
\value{
	\item{design}{(Partially balanced semi-Latin rectangle of cell size two for v treatments}
	\item{Avg.Effi}{Average efficiency of the design}
  }
\examples{
pbslr(5)
pbslr(6)
}
\author{Rajender Parsad <rajender.parsad@icar.gov.in>}
\keyword{partially balanced}
\keyword{semi-Latin rectangle}