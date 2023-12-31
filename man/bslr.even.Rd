\name{bslr.even}
\alias{bslr.even}
\title{Balanced semi-Latin rectangles for given cell size and even number of treatments}
\description{This function generates a balanced semi-Latin rectangles for given cell size and even number of treatments}
\usage{bslr.even(v, k)}
\arguments{
	\item{v}{number of treatments and is even}
	\item{k}{cell size, preferably up to 10}
	}
\value{
	\item{design}{(Balanced semi-Latin rectangle of cell size two for v treatments}
	\item{Avg.Effi}{Average efficiency of the design}
  }
\examples{
bslr.even(4,2)
}
\author{Sukanta Dash <sukanta.dash@icar.gov.in>}
\keyword{balanced}
\keyword{semi-Latin rectangle}