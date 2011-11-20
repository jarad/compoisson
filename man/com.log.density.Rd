\name{com.log.density}
\alias{com.log.density}
\title{ Computes the Log PMF of the COM-Poisson Distribution }
\description{
	Computes the log probability mass function of the COM-Poisson distribution for
	given values of the parameters.
}
\usage{
com.log.density(x, lambda, nu, log.z = NULL)
}
\arguments{
	\item{x}{ level to evaulate the log PMF at }
	\item{lambda}{ value of the lambda parameter }
	\item{nu}{ value of the nu parameter }
	\item{log.z}{ log of the normalizing constant, computed if not specified }
}
\details{
	Computes the log probability mass function of the COM-Poisson distribution
	\deqn{
		\log f(x) = x \log \lambda - \log(Z(\lambda,\nu)) - \nu \sum_{i=1}^x x.
	}{log f(x) = x * log(lambda) - log(Z) - nu * log(x!). }
}
\value{
	The log probability that a random COM-Poisson variable X takes value x.
}
\references{ Shmueli, G., Minka, T. P., Kadane, J. B., Borle, S. and Boatwright, P., \dQuote{A useful distribution for fitting discrete data: Revival of the Conway-Maxwell-Poisson distribution,} J. Royal Statist. Soc., v54, pp. 127-142, 2005. }
\author{ Jeffrey Dunn }
\seealso{ \code{\link{com.loglikelihood}}, \code{\link{dcom}} }
\examples{
	data(insurance);
	fit = com.fit(Lemaire);
	com.log.density(0, fit$lambda, fit$nu, fit$z);
}
\keyword{ models }
