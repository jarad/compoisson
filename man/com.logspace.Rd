\name{com.logspace}
\alias{com.log.difference}
\alias{com.log.sum}
\alias{com.log.factorial}
\title{ Operations in Log-space }
\description{
	Computes the difference of two values in log-space.
}
\usage{
	com.log.difference(x, y)
	com.log.sum(x, y)
	com.log.factorial(x)
}
\arguments{
  \item{x}{ first value }
  \item{y}{ second value }
}
\details{
	\code{com.log.difference} computes the difference of two values in log-space, \eqn{log( e^x - e^y )},
	without significant	chance of overflow or underflow.

	\code{com.log.sum} computes the sum of two values in log-space, \eqn{log( e^x + e^y )}, without
	significant change of overflow or underflow.

	\code{com.log.factorial} computes \eqn{log(x!)} which is equivalent to a summation.
}
\value{
	The requested computation in log-space.
}
\author{ Jeffrey Dunn }
\examples{
	a = exp(com.log.difference(log(100), log(20))); # a = 80
	b = exp(com.log.sum(log(100), log(20))); # b = 120
	c = exp(com.log.factorial(4)); # c = 24
}
\keyword{ manip }
