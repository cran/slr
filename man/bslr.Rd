\name{bslr}
\alias{bslr}
\title{Balanced semi-Latin rectangles for given number of treatments and cell size}
\description{This function generates a balanced semi-Latin rectangles for given number of treatments and cell size}
\usage{bslr(v, k)}
\arguments{
	\item{v}{number of treatments}	
	\item{k}{cell size, preferably up to 10}
	}
\value{
	\item{design}{(Balanced semi-Latin rectangle of cell size k for v treatments}
	\item{Avg.Effi}{Average efficiency of the design}
  }
\examples{
bslr(5, 2)
bslr(7, 3)
}
\author{Kaushal Kumar Yadav <kaushalyadav0796@gmail.com>}
\keyword{balanced}
\keyword{semi-Latin rectangle}