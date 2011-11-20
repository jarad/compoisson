\name{com.confint}
\alias{com.confint}
\title{ Computes a confidence interval for parameter estimates of the COM-Poisson Distribution }
\description{
	Computes a pivotal bootstrap confidence interval for maximum likelihood parameter
	estimates.
}
\usage{
	com.confint(data, level=0.95, B=1000, n=1000)
}
\arguments{
	\item{data}{ the matrix of data to fit }
	\item{level}{ the level of the confidence interval }
	\item{B}{ number of repetitions of the bootstrap }
	\item{n}{ number of data points in each bootstrap sample }
}
\details{
	Uses a standard pivotal confidence interval from a bootstrap sample.
}
\value{
	A matrix containing the confidence intervals for each parameter
}
\references{ Wasserman, L. (2005). "All of Statistics: A Concise Course in Statistical Inference," Springer Texts in Statistics. }
\author{ Akshaya Jha, Jeffrey Dunn }
\seealso{ \code{\link{com.fit}} }
\keyword{ models }
