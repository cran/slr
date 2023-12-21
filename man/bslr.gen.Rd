\name{bslr.gen}
\alias{bslr.gen}
\title{Balanced semi-Latin rectangles for given number of treatments and cell size}
\description{This function generates a balanced semi-Latin rectangles for given number of treatments and cell size}
\usage{bslr.gen(v, k)}
\arguments{
	\item{v}{number of treatments}	
	\item{k}{cell size, preferably up to 10}
	}
\value{
	\item{design}{(Balanced semi-Latin rectangle of cell size k with v treatments}
	\item{Avg.Effi}{Average efficiency of the design}
  }
\examples{
bslr.gen(7, 4)
}
\author{Kaushal Kumar Yadav <kaushalyadav0796@gmail.com>}
\keyword{balanced}
\keyword{semi-Latin rectangle}