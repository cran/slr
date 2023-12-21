\name{pbslr}
\alias{pbslr}
\title{Partially balanced semi-Latin rectangles of cell size two, three and four}
\description{This function generates a partially balanced semi-Latin rectangles of cell size two, three and four}
\usage{pbslr(v,k)}
\arguments{
	\item{v}{number of treatments}	
	\item{k}{cell size, 2, 3 or 4 is supported.}
	}
\value{
	\item{design}{(Partially balanced semi-Latin rectangle of cell size two, three or four for v treatments}
	\item{Avg.Effi}{Average efficiency of the design}
  }
\examples{
pbslr(5, 2)
pbslr(6, 3)
}
\author{Rajender Parsad <rajender.parsad@icar.gov.in>}
\keyword{partially balanced}
\keyword{semi-Latin rectangle}