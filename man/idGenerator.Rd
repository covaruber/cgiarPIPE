\name{idGenerator}
\alias{idGenerator}

\title{
ID generator
}
\description{
Generates an id to be assigned to each analysis.
}
\usage{
idGenerator(nn=5,nl=5)
}

\arguments{
  
\item{nn}{
number of digits for the code.
}
\item{nl}{
number of letters for the code.
}

}
\details{

At high marker density, the relationship matrix is estimated as \eqn{A=W W'/c}, where \eqn{W_{ik} = X_{ik} + 1 - 2 p_k} and \eqn{p_k} is the frequency of the 1 allele at marker k.  By using a normalization constant of \eqn{c = 2 \sum_k {p_k (1-p_k)}}, the mean of the diagonal elements is \eqn{1 + f} (Endelman and Jannink 2012).

}
\value{
If return.imputed = FALSE, the \eqn{n \times n} additive relationship matrix is returned.

If return.imputed = TRUE, the function returns a list containing
\describe{
\item{$out}{A character vector of same length as \code{Source} in which the \code{Search}-values have been replaced by the corresponding \code{Replace}-values.}
}
}

\references{
None
}

\examples{

idGenerator()

}