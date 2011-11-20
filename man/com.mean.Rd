\name{com.mean}
\alias{com.mean}
\title{ Computes Mean of the COM-Poisson Distribution  }
\description{
	Computes the mean of the COM-Poisson distribution for given values
	of the parameters.
}
\usage{
	com.mean(lambda, nu)
}
\arguments{
	\item{lambda}{ value of lambda parameter }
	\item{nu}{ value of the nu parameter }
}
\details{
	Uses \code{\link{com.expectation}} to compute the first moment of the distribution.
}
\value{
	The mean of the distribution.
}
\references{ Shmueli, G., Minka, T. P., Kadane, J. B., Borle, S. and Boatwright, P., \dQuote{A useful distribution for fitting discrete data: Revival of the Conway-Maxwell-Poisson distribution,} J. Royal Statist. Soc., v54, pp. 127-142, 2005. }
\author{ Jeffrey Dunn }
\seealso{ \code{\link{com.expectation}}, \code{\link{com.var}} }
\examples{
	data(insurance)
	model = com.fit(Lemaire);
	com.mean(model$lambda, model$nu);
}
\keyword{ models }
