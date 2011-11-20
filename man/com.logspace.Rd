\name{com.logspace}
\alias{com.log.sum}
\title{ Operations in Log-space }
\description{
	Computes the difference of two values in log-space.
}
\usage{
	com.log.sum(x, y)
}
\arguments{
  \item{x}{ first value }
  \item{y}{ second value }
}
\details{
	\code{com.log.sum} computes the sum of two values in log-space, \eqn{log( e^x + e^y )}, without
	significant change of overflow or underflow.
}
\value{
	The requested computation in log-space.
}
\author{ Jeffrey Dunn }
\examples{
	b = exp(com.log.sum(log(100), log(20))); # b = 120
}
\keyword{ manip }
