\name{compoisson-data}
\alias{insurance}
\alias{Lemaire}
\alias{Buhlmann}
\docType{data}
\title{ Insurance Count Datasets }
\description{
	Two auto insurance datasets compiled from published works. The Lemaire dataset
	contains published aggregate claim numbers for automobile third-party liability
	insurance of a Belgian insurance company in the early 1990's. The Buhlmann dataset
	originates from aggregate accident claims in 1961 for a class of auto insurance in
	Switzerland.
}
\usage{data(insurance)}
\format{
	Each dataset is a matrix with two columns. The first column contains the
	levels and the second contains the number of customers who submitted the
	corresponding level of claims.
}
\source{
	Lemaire, Jean. \dQuote{Bonus-Malus Systems for Automobile Insurance}. Kluwer Academic Publishers, 1995.

	Panjer, Harry. \dQuote{Actuarial Mathematics (Proceedings of Symposia in Applied Mathematics)}. Providence: American Mathematical Society, 1986.
	
}
\examples{
	data(insurance)
	Lemaire
	Buhlmann
}
\keyword{datasets}
