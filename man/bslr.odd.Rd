\name{bslr.odd}
\alias{bslr.odd}
\title{Balanced semi-Latin rectangles of for given cell size and odd number of treatments}
\description{This function generates a partially balanced semi-Latin rectangles for given cell size and odd number of treatments}
\usage{bslr.odd(v,k)}
\arguments{
	\item{v}{number of treatments and is odd}
	\item{k}{cell size, prerably up to 10}
	}
\value{
	\item{design}{(Balanced semi-Latin rectangle of cell size two for v treatments}
	\item{Avg.Effi}{Average efficiency of the design}
  }
\examples{
bslr.odd(5,2)
}
\author{Kaushal Kumar Yadav <kaushalyadav0796@gmail.com>}
\keyword{balanced}
\keyword{semi-Latin rectangle}