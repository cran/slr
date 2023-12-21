\name{pbslr.even}
\alias{pbslr.even}
\title{Partially balanced semi-Latin rectangles of cell size two, three and four for even number of teratments}
\description{This function generates a partially balanced semi-Latin rectangles of cell size two, three and four for even number of treatments}
\usage{pbslr.even(v,k)}
\arguments{
	\item{v}{number of treatments and is even}
	\item{k}{cell size, 2, 3 and 4 is supported}	
	}
\value{
	\item{design}{(Partially balanced semi-Latin rectangle of cell size two, three and four for v treatments for even number of treatments}
	\item{Avg.Effi}{Average efficiency of the design}
  }
\examples{
pbslr.even(4,2)
pbslr.even(6,4)
}
\author{Baidya Nath Mandal <mandal.stat@gmail.com>}
\keyword{partially balanced}
\keyword{semi-Latin rectangle}