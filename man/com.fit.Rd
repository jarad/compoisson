\name{com.fit}
\alias{com.fit}
\title{ Computes COM-Poisson Regression }
\description{
	Computes the maximum likelihood estimates of the COM-Poisson model for given count data.
}
\usage{com.fit(x)}
\arguments{
	\item{x}{ matrix of count data }
}
\details{
	The argument x should consist of a matrix where the first column is the level and the second column is the count for the corresponding level.
}
\value{
	Returns an object containing four fields:
 	\item{lambda }{ Estimate of the lambda parameter}
	\item{nu }{ Estimate of the nu parameter}
	\item{z }{ Normalizing constant}
	\item{fitted.values }{ Estimated counts at given levels }
}
\references{ Shmueli, G., Minka, T. P., Kadane, J. B., Borle, S. and Boatwright, P., \dQuote{A useful distribution for fitting discrete data: Revival of the Conway-Maxwell-Poisson distribution,} J. Royal Statist. Soc., v54, pp. 127-142, 2005. }
\author{ Jeffrey Dunn }
\seealso{ \code{\link{com.compute.z}}, \code{\link{com.loglikelihood}} }
\examples{
	data(insurance)
	com.fit(Lemaire);
}

\keyword{ models }
\keyword{ regression }
