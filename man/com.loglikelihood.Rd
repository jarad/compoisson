\name{com.loglikelihood}
\alias{com.loglikelihood}
\title{ Computes Log-Likelihood of COM-Poisson }
\description{
	Given a set of data, computes the log-likelihood of the data under the
	COM-Poisson distribution for values of the parameters.
}
\usage{
com.loglikelihood(x, lambda, nu)
}
\arguments{
	\item{x}{ matrix of count data }
	\item{lambda}{ value of lambda parameter }
	\item{nu}{ value of nu parameter }
}
\details{
	The argument x should consist of a matrix where the first column is the level and the second column is the count for the corresponding level.
}
\value{
	The log-likelihood of the data.
}
\references{ Shmueli, G., Minka, T. P., Kadane, J. B., Borle, S. and Boatwright, P., \dQuote{A useful distribution for fitting discrete data: Revival of the Conway-Maxwell-Poisson distribution,} J. Royal Statist. Soc., v54, pp. 127-142, 2005. }
\author{ Jeffrey Dunn }
\seealso{ \code{\link{com.fit}}, \code{\link{dcom}} }

\keyword{ models }
