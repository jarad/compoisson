\name{dcom}
\alias{dcom}
\alias{rcom}
\alias{qcom}
\title{ The COM-Poisson Distribution }
\description{
	Probability mass function and random generation for the COM-Poisson
	distribution for given values of the parameters.
}
\usage{
	dcom(x, lambda, nu, z = NULL)
        qcom(p, lambda, nu, log.z = NULL, lower.tail = TRUE, log.p = FALSE)
	rcom(n, lambda, nu, log.z = NULL)
}
\arguments{
	\item{x}{ vector of (non-negative integer) quantiles. }
        \item{p}{ vector of probabilities.}
	\item{lambda}{ scalar value of lambda parameter }
	\item{nu}{ scalar value of nu parameter }
	\item{z}{ normalizing constant, computed if not specified }
	\item{n}{ number of random values to return }
	\item{log.z}{ natural log of z }
        \item{log.p}{ if TRUE, probabilities p are given as log(p).}
        \item{lower.tail}{ if TRUE (default), probabilities are P[X <= x],
          otherwise, P[X > x].}
}
\details{
	Computes the probability mass function of the COM-Poisson distribution
	\deqn{
		f(x) = \frac{1}{Z(\lambda,\nu)}\frac{\lambda^x}{(x!)^\nu}.
	}{f(x) = (1/Z) (lambda^x)/(x!^nu).}
}
\value{
	\code{dcom} gives the probability that a random COM-Poisson variable X takes value x.
	\code{rcom} gives a vector of \code{n} random values sampled from the COM-Poisson
	distribution.
        \code{qcom} gives the quantile function
}
\references{ Shmueli, G., Minka, T. P., Kadane, J. B., Borle, S. and Boatwright, P., \dQuote{A useful distribution for fitting discrete data: Revival of the Conway-Maxwell-Poisson distribution,} J. Royal Statist. Soc., v54, pp. 127-142, 2005. }
\author{ Jeffrey Dunn }
\seealso{ \code{\link{com.loglikelihood}}, \code{\link{com.log.density}} }
\examples{
	data(insurance);
	fit = com.fit(Lemaire);
	dcom(0, fit$lambda, fit$nu, fit$z);
	r = rcom(10, fit$lambda, fit$nu);
}

\keyword{ models }
