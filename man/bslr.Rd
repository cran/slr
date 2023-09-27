\name{bslr}
\alias{bslr}
\title{Balanced semi-Latin rectangles of cell size two and three}
\description{This function generates a balanced semi-Latin rectangles of cell size two and three}
\usage{bslr(v, k)}
\arguments{
	\item{v}{number of treatments}	
	\item{k}{cell size. Either 2 or 3 is supported.}
	}
\value{
	\item{design}{(Balanced semi-Latin rectangle of cell size two for v treatments}
	\item{Avg.Effi}{Average efficiency of the design}
  }
\examples{
bslr(5, 2)
bslr(7, 3)
}
\author{Kaushal Kumar Yadav <kaushalyadav0796@gmail.com>}
\keyword{balanced}
\keyword{semi-Latin rectangle}