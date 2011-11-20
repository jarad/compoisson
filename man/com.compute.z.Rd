\name{com.compute.z}
\alias{com.compute.z}
\alias{com.compute.log.z}
\title{ Compute COM-Poisson Normalizing Constant }
\description{
	Computes the normalizing constant in the COM-Poisson model for given values
	of the parameters.
}
\usage{
	com.compute.z(lambda, nu, log.error = 0.001)
	com.compute.log.z(lambda, nu, log.error = 0.001)
}
\arguments{
	\item{lambda}{ Lambda value in COM-Poisson distribution }
	\item{nu}{ Nu value in COM-Poisson distribution }
	\item{log.error}{ Precision in the log of the normalizing constant }
}
\details{
	\code{com.compute.z} computes the COM-Poisson normalizing constant
	\deqn{
		z = \sum_{i=0}^\infty \frac{\lambda^j}{(j!)^\nu}
	}{z = Sum (lambda^j)/(j!^nu) }
	to the specified precision. If no precision is specified, then the package default
	is used.

	\code{com.compute.log.z} is equivalent to \code{log(com.compute.z(lambda, nu))} but
	provudes additional precision.
}
\value{
	The normalizing constant as a real number with specified precision.
}
\references{ Shmueli, G., Minka, T. P., Kadane, J. B., Borle, S. and Boatwright, P., \dQuote{A useful distribution for fitting discrete data: Revival of the Conway-Maxwell-Poisson distribution,} J. Royal Statist. Soc., v54, pp. 127-142, 2005. }
\author{ Jeffrey Dunn }
\seealso{ \code{\link{com.fit}} }
\examples{
	data(insurance);
	fit = com.fit(Lemaire);
	z = com.compute.z(fit$lambda, fit$nu);
}
\keyword{ models }
