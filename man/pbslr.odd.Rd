\name{pbslr.odd}
\alias{pbslr.odd}
\title{Partially balanced semi-Latin rectangles of cell size two and three for odd number of teratments}
\description{This function generates a partially balanced semi-Latin rectangles of cell size two and three for odd number of treatments}
\usage{pbslr.odd(v,k)}
\arguments{
	\item{v}{number of treatments and is odd}
	\item{k}{cell size. Either 2 or 3 is supported}	
	}
\value{
	\item{design}{(Partially balanced semi-Latin rectangle of cell size two for v treatments for odd number of treatments}
	\item{Avg.Effi}{Average efficiency of the design}
  }
\examples{
pbslr.odd(5,2)
}
\author{Baidya Nath Mandal <mandal.stat@gmail.com>}
\keyword{partially balanced}
\keyword{semi-Latin rectangle}