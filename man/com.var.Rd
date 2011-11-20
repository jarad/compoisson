\name{com.var}
\alias{com.var}
\title{ Computes Variance of the COM-Poisson Distribution  }
\description{
	Computes the variance of the COM-Poisson distribution for given values
	of the parameters.
}
\usage{
	com.var(lambda, nu)
}
\arguments{
	\item{lambda}{ value of lambda parameter }
	\item{nu}{ value of the nu parameter }
}
\details{
	Uses \code{\link{com.expectation}} to compute the second moment of the distribution
	and subtracts the squared mean, computed using \code{\link{com.mean}}.
}
\value{
	The variance of the distribution.
}
\references{ Shmueli, G., Minka, T. P., Kadane, J. B., Borle, S. and Boatwright, P., \dQuote{A useful distribution for fitting discrete data: Revival of the Conway-Maxwell-Poisson distribution,} J. Royal Statist. Soc., v54, pp. 127-142, 2005. }
\author{ Jeffrey Dunn }
\seealso{ \code{\link{com.expectation}}, \code{\link{com.mean}} }
\examples{
	data(insurance)
	model = com.fit(Lemaire);
	com.var(model$lambda, model$nu);
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
