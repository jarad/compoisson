\name{com.expectation}
\alias{com.expectation}
\title{ Computes Expectation of a Function of a COM-Poisson Random Variable }
\description{
	Computes an expectation of a function of a COM-Poisson random variable.
}
\usage{
com.expectation(f, lambda, nu, log.error = 0.001)
}
\arguments{
	\item{f}{ function taking as a single argument the value of x }
	\item{lambda}{ value of lambda parameter }
	\item{nu}{ value of nu parameter }
	\item{log.error}{ precision in the log of the expectation }
}
\details{
	Computes the expectation \eqn{E[f(X)]}{E[f(X)]} where X is a COM-Poisson random variable.
}
\value{
	The expectation as a real number.
}
\references{ Shmueli, G., Minka, T. P., Kadane, J. B., Borle, S. and Boatwright, P., \dQuote{A useful distribution for fitting discrete data: Revival of the Conway-Maxwell-Poisson distribution,} J. Royal Statist. Soc., v54, pp. 127-142, 2005. }
\author{ Jeffrey Dunn }
\seealso{ \code{\link{com.mean}}, \code{\link{com.var}}, \code{\link{com.fit}} }
\keyword{ models }
